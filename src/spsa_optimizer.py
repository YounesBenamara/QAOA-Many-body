"""
Core SPSA Optimizer and utility functions for VQE on hardware/simulators.

"""

import time
import json
import datetime as dt
from pathlib import Path
import numpy as np

import sys
import platform
from importlib.metadata import version as pkg_version

# ═══════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════

def utc_now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()

# Extracts a standard Python float from a numpy array or Qiskit EV output
def scalar(x) -> float:
    return float(np.asarray(x).reshape(-1)[0])

def wrap_angles(theta: np.ndarray) -> np.ndarray:
    return (theta + np.pi) % (2 * np.pi) - np.pi

def window_mean_and_se(values, stds):
    values = np.asarray(values, dtype=float)
    stds = np.asarray(stds, dtype=float)
    mu = float(np.mean(values))
    se = float(np.sqrt(np.sum(stds**2)) / len(stds))
    return mu, se

def plateau_test(energies, stds, window, z, deterministic_tol=1e-4):
    """
    Noisy mode:
        stop if improvement < z * combined_SE

    Deterministic fallback:
        if threshold == 0, stop if improvement < deterministic_tol
    """
    if len(energies) < 2 * window:
        return False, {}

    prev_vals = np.asarray(energies[-2 * window:-window], dtype=float)
    curr_vals = np.asarray(energies[-window:], dtype=float)
    prev_stds = np.asarray(stds[-2 * window:-window], dtype=float)
    curr_stds = np.asarray(stds[-window:], dtype=float)

    mu_prev, se_prev = window_mean_and_se(prev_vals, prev_stds)
    mu_curr, se_curr = window_mean_and_se(curr_vals, curr_stds)

    improvement = mu_prev - mu_curr
    threshold = z * np.sqrt(se_prev**2 + se_curr**2)

    if threshold == 0.0:
        stop = improvement < deterministic_tol
    else:
        stop = improvement < threshold

    info = {
        "mu_prev": mu_prev,
        "mu_curr": mu_curr,
        "se_prev": se_prev,
        "se_curr": se_curr,
        "improvement": improvement,
        "threshold": threshold,
    }
    return stop, info

def to_jsonable(obj):
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(x) for x in obj]
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if hasattr(obj, "isoformat"):
        try:
            return obj.isoformat()
        except Exception:
            pass
    return str(obj)

def append_jsonl(path: Path, record: dict) -> None:
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(to_jsonable(record), ensure_ascii=False) + "\n")

def save_json(path: Path, data: dict) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(to_jsonable(data), f, ensure_ascii=False, indent=2)


#---ESTIMATOR TRIPLET---

# Submits exactly three executions to the Estimator in a single batch to calculate the gradient and energy
def estimate_triplet(estimator, circuit, observable, theta_plus, theta_minus, theta_center):
    pubs = [
        (circuit, observable, theta_plus),
        (circuit, observable, theta_minus),
        (circuit, observable, theta_center),
    ]
    t0 = time.time()
    job = estimator.run(pubs)
    result = job.result()
    t1 = time.time()

    pub_plus = result[0]
    pub_minus = result[1]
    pub_center = result[2]

    e_plus = scalar(pub_plus.data.evs)
    s_plus = scalar(pub_plus.data.stds)

    e_minus = scalar(pub_minus.data.evs)
    s_minus = scalar(pub_minus.data.stds)

    e_center = scalar(pub_center.data.evs)
    s_center = scalar(pub_center.data.stds)
    
    # Handle the difference in job_id logic between real IBM connection and local FakeBackend
    try:
        if callable(getattr(job, "job_id", None)):
            job_id = job.job_id()
        else:
            job_id = getattr(job, "job_id", "fake_local")       
    except Exception:
        job_id = "fake_local"

    meta = {
        "job_id": job_id,
        "wall_time_s": t1 - t0,
        "pub_plus_metadata": getattr(pub_plus, "metadata", {}),
        "pub_minus_metadata": getattr(pub_minus, "metadata", {}),
        "pub_center_metadata": getattr(pub_center, "metadata", {}),
    }

    return (e_plus, s_plus), (e_minus, s_minus), (e_center, s_center), meta


#---INITIAL PARAMS---

# Determines starting parameters: resumes from snapshot if available, or adds small random noise around zero
def choose_initial_params(n_params: int, snapshot_npz: Path | None, init_mode: str, init_std_dev: float, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    if init_mode == "resume" and snapshot_npz is not None and snapshot_npz.exists():
        data = np.load(snapshot_npz)
        theta = np.asarray(data["theta_opt"], dtype=float)
        print(f"Loaded warm-start parameters from {snapshot_npz}")
        return theta
        

    # We might change the standard dev for a warm start
    theta = init_std_dev * rng.standard_normal(n_params)
    print("Starting from random initial parameters")
    return theta

    
#---VQE DRIVER---

def run_vqe_spsa_decay_hardware(ansatz_isa,hamiltonian_isa,estimator,
                                initial_params,run_dir: Path,
                                maxiter=30, a_lr=0.08,a_shift=10.0,
                                alpha=0.602,c_pert=0.10,gamma=0.101,
                                plateau_window=10, plateau_z=1.5,
                                plateau_patience=3,tail_window=37,
                                enable_early_stop=True, resume_history=False,
                                seed=42,
):
    """
    Core manual SPSA loop handling decaying gain sequences, 
    early stopping plateaus, and incremental data logging.

    Energy reporting uses Polyak-Ruppert tail averaging (mean of last
    `tail_window` iterations) for a statistically robust estimate.
    Parameter selection keeps the running-best theta (lowest window mean).
    """
    rng_local = np.random.default_rng(seed)

    history_jsonl = run_dir / "history.jsonl"
    snapshot_npz = run_dir / "latest_snapshot.npz"

    theta = np.asarray(initial_params, dtype=float).copy()
    theta = wrap_angles(theta)

    n_params = theta.size

    # --- Resume: load previous histories to concatenate ---
    prior_iters = 0
    if resume_history and snapshot_npz.exists():
        prev = np.load(snapshot_npz, allow_pickle=True)
        cost_history    = list(prev["cost_history"])
        std_history     = list(prev["std_history"])
        grad_norm_history = list(prev["grad_norm_history"])
        ak_history      = list(prev["ak_history"])
        ck_history      = list(prev["ck_history"])
        job_ids         = list(prev["job_ids"])
        prior_iters     = len(cost_history)
        print(f"[RESUME] Loaded {prior_iters} previous iterations from snapshot")
    else:
        cost_history = []
        std_history = []
        grad_norm_history = []
        ak_history = []
        ck_history = []
        job_ids = []

    plateau_hits = 0
    converged_early = False

    best_window_mean = np.inf
    best_window_se = np.inf
    best_theta = theta.copy()

    for k in range(maxiter):
        k_global = k + prior_iters  # continue decay from previous run
        delta = rng_local.choice([-1.0, 1.0], size=n_params)

        ak = a_lr / ((k_global + 1 + a_shift) ** alpha)
        ck = c_pert / ((k_global + 1) ** gamma)

        theta_center = theta.copy()
        theta_plus = wrap_angles(theta_center + ck * delta)
        theta_minus = wrap_angles(theta_center - ck * delta)

        (e_plus, s_plus), (e_minus, s_minus), (e_center, s_center), meta = estimate_triplet(
            estimator,
            ansatz_isa,
            hamiltonian_isa,
            theta_plus,
            theta_minus,
            theta_center,
        )

        grad = ((e_plus - e_minus) / (2.0 * ck)) * delta
        grad_norm = float(np.linalg.norm(grad))

        cost_history.append(e_center)
        std_history.append(s_center)
        grad_norm_history.append(grad_norm)
        ak_history.append(ak)
        ck_history.append(ck)
        job_ids.append(meta.get("job_id", "fake_local"))

        if len(cost_history) >= plateau_window:
            mu_curr, se_curr = window_mean_and_se(
                cost_history[-plateau_window:],
                std_history[-plateau_window:],
            )
            if mu_curr < best_window_mean:
                best_window_mean = mu_curr
                best_window_se = se_curr
                best_theta = theta_center.copy()

        stop_now, info = plateau_test(
            cost_history,
            std_history,
            window=plateau_window,
            z=plateau_z,
        )

        if stop_now:
            plateau_hits += 1
        else:
            plateau_hits = 0

        # during the run print the current energy and standard deviation
        msg = f"Iter {k_global+1:03d}/{prior_iters+maxiter} | E={e_center:+.6f} ± {s_center:.6f}"
        if plateau_hits > 0 and info:
            msg += (
                f"| Δwindow={info['improvement']:+.6f} vs thr={info['threshold']:.6f} "
                f"| hits={plateau_hits}/{plateau_patience}"
            )
        print(msg)

        # Incremental save
        iteration_record = {
            "iteration": k_global + 1,
            "timestamp_utc": utc_now_iso(),
            "job_id": meta.get("job_id", "fake_local"),
            "job_wall_time_s": meta["wall_time_s"],
            "theta_center": theta_center,
            "theta_plus": theta_plus,
            "theta_minus": theta_minus,
            "delta": delta,
            "ak": ak,
            "ck": ck,
            "e_plus": e_plus,
            "std_plus": s_plus,
            "e_minus": e_minus,
            "std_minus": s_minus,
            "e_center": e_center,
            "std_center": s_center,
            "grad_norm": grad_norm,
            "plateau_hits": plateau_hits,
            "plateau_info": info,
            "pub_plus_metadata": meta["pub_plus_metadata"],
            "pub_minus_metadata": meta["pub_minus_metadata"],
            "pub_center_metadata": meta["pub_center_metadata"],
        }
        append_jsonl(history_jsonl, iteration_record)

        np.savez(
            snapshot_npz,
            theta_current=theta_center,
            theta_opt=best_theta,
            energy_opt=best_window_mean if np.isfinite(best_window_mean) else e_center,
            energy_opt_se=best_window_se if np.isfinite(best_window_se) else s_center,
            cost_history=np.asarray(cost_history, dtype=float),
            std_history=np.asarray(std_history, dtype=float),
            grad_norm_history=np.asarray(grad_norm_history, dtype=float),
            ak_history=np.asarray(ak_history, dtype=float),
            ck_history=np.asarray(ck_history, dtype=float),
            job_ids=np.asarray(job_ids, dtype=object),
            converged_early=converged_early,
            last_iteration=k_global + 1,
        )

        if enable_early_stop and plateau_hits >= plateau_patience:
            converged_early = True
            break

        theta = wrap_angles(theta_center - ak * grad)

    # --- Polyak-Ruppert tail averaging for energy estimate ---
    n_done = len(cost_history)
    tw = min(tail_window, n_done)  # guard against short runs
    pr_energy, pr_se = window_mean_and_se(
        cost_history[-tw:],
        std_history[-tw:],
    )

    # --- Parameters: keep running-best theta ---
    if not np.isfinite(best_window_mean):
        best_theta = theta_center.copy()

    return {
        "theta_opt": best_theta,
        "energy_opt": pr_energy,
        "energy_opt_se": pr_se,
        "energy_best_window": best_window_mean if np.isfinite(best_window_mean) else cost_history[-1],
        "energy_best_window_se": best_window_se if np.isfinite(best_window_se) else std_history[-1],
        "cost_history": np.asarray(cost_history, dtype=float),
        "std_history": np.asarray(std_history, dtype=float),
        "grad_norm_history": np.asarray(grad_norm_history, dtype=float),
        "ak_history": np.asarray(ak_history, dtype=float),
        "ck_history": np.asarray(ck_history, dtype=float),
        "job_ids": job_ids,
        "converged_early": converged_early,
        "tail_window_used": tw,
    }


#---DATA LOGGING UTILS---

def build_and_save_metadata(
    run_dir, backend, auth_mode,
    N, g, N_LAYERS, REAL_ANSATZ,
    SPSA_MAXITER, A_LR, A_SHIFT, ALPHA, C_PERT, GAMMA,
    SHOTS, USE_RESILIENCE_LEVEL_1,
    PLATEAU_WINDOW, PLATEAU_Z, PLATEAU_PATIENCE,
    INIT_MODE, INIT_STD_DEV, initial_params,
    E_exact, ansatz, ansatz_isa
):
    metadata = {
        "timestamp_utc": utc_now_iso(),
        "auth_mode": auth_mode,
        "backend_name": backend.name,
        "backend_num_qubits": getattr(backend, "num_qubits", None),
        "backend_version": str(getattr(backend, "backend_version", "")),
        "physics": {
            "N": N, "g": g, "n_layers": N_LAYERS, "real_ansatz": REAL_ANSATZ,
        },
        "optimizer": {
            "name": "spsa_decay_manual",
            "maxiter": SPSA_MAXITER, "A_LR": A_LR, "A_SHIFT": A_SHIFT, "ALPHA": ALPHA, "C_PERT": C_PERT, "GAMMA": GAMMA,
        },
        "runtime": {
            "shots": SHOTS, "use_resilience_level_1": USE_RESILIENCE_LEVEL_1, "job_mode": True, "session_mode": False,
        },
        "plateau": {
            "window": PLATEAU_WINDOW, "z": PLATEAU_Z, "patience": PLATEAU_PATIENCE,
        },
        "initialization": {
            "mode": INIT_MODE, "std_dev": INIT_STD_DEV, "initial_params": initial_params,
        },
        "reference": {
            "exact_ground_energy": E_exact,
        },
        "circuits": {
            "logical_depth": ansatz.depth(),
            "logical_gate_count": dict(ansatz.count_ops()),
            "isa_depth": ansatz_isa.depth(),
            "isa_gate_count": dict(ansatz_isa.count_ops()),
            "layout": str(ansatz_isa.layout),
            "final_index_layout": to_jsonable(
                ansatz_isa.layout.final_index_layout(filter_ancillas=True)
            ) if hasattr(ansatz_isa.layout, "final_index_layout") else None,
        },
        "software": {
            "python": sys.version, "platform": platform.platform(), "qiskit": pkg_version("qiskit"), "qiskit_ibm_runtime": pkg_version("qiskit-ibm-runtime"), "numpy": pkg_version("numpy"),
        },
    }
    save_json(run_dir / "metadata.json", metadata)


def save_summary_and_arrays(
    run_dir, backend_name, result, E_exact, wall_time_total_s
):
    theta_opt = result["theta_opt"]
    E_vqe = float(result["energy_opt"])
    E_vqe_se = float(result["energy_opt_se"])
    cost_history = result["cost_history"]
    std_history = result["std_history"]
    grad_norm_history = result["grad_norm_history"]
    ak_history = result["ak_history"]
    ck_history = result["ck_history"]
    job_ids = result["job_ids"]
    converged_early = bool(result["converged_early"])
    
    rel_error = abs(E_vqe - E_exact) / abs(E_exact) if E_exact != 0 else 0

    summary = {
        "timestamp_utc": utc_now_iso(),
        "run_dir": str(run_dir),
        "backend_name": backend_name,
        "iterations_done": len(cost_history),
        "converged_early": converged_early,
        "exact_energy": E_exact,
        "vqe_energy": E_vqe,
        "vqe_energy_se": E_vqe_se,
        "relative_error": rel_error,
        "job_ids": job_ids,
        "wall_time_total_s": wall_time_total_s,
        "theta_opt": theta_opt,
        "cost_history": cost_history,
        "std_history": std_history,
        "grad_norm_history": grad_norm_history,
        "ak_history": ak_history,
        "ck_history": ck_history,
    }
    save_json(run_dir / "summary.json", summary)

    np.savez(
        run_dir / "final_arrays.npz",
        theta_opt=np.asarray(theta_opt, dtype=float),
        cost_history=np.asarray(cost_history, dtype=float),
        std_history=np.asarray(std_history, dtype=float),
        grad_norm_history=np.asarray(grad_norm_history, dtype=float),
        ak_history=np.asarray(ak_history, dtype=float),
        ck_history=np.asarray(ck_history, dtype=float),
        exact_energy=np.asarray([E_exact], dtype=float),
        final_energy=np.asarray([E_vqe], dtype=float),
        final_energy_se=np.asarray([E_vqe_se], dtype=float),
        converged_early=np.asarray([converged_early]),
    )
    
    return rel_error, len(cost_history), converged_early


def print_console_summary(backend_name, E_exact, E_vqe, E_vqe_se, rel_error, iterations, converged_early, wall_time_s, run_dir):
    print("\n" + "─" * 60)
    print(f"Backend         : {backend_name}")
    print(f"Exact energy    : {E_exact:.6f}")
    print(f"VQE energy      : {E_vqe:.6f} ± {E_vqe_se:.6f}")
    print(f"Relative error  : {rel_error*100:.3f} %")
    print(f"Iterations      : {iterations}")
    print(f"Early stop      : {converged_early}")
    print(f"Wall time total : {wall_time_s:.1f} s")
    print(f"Saved in        : {run_dir}")
    print("─" * 60)
