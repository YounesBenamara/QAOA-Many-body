"""
VQE on FakeMarrakesh — mirrors the real-hardware pipeline (vqe_hardware.py).

Purpose
-------
- Run the exact same VQE workflow as vqe_hardware.py but on FakeMarrakesh
  (noisy classical simulator that mimics ibm_marrakesh).
- Same transpilation, same Estimator(mode=backend), same SPSA with decay,
  same incremental saving (JSONL + NPZ + metadata.json + summary.json).
- This lets you validate the full pipeline locally before spending QPU time.
"""

from __future__ import annotations

import time
import datetime as dt
from pathlib import Path

from qiskit import qpy
from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit_ibm_runtime.fake_provider import FakeMarrakesh
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit.visualization import plot_gate_map
import matplotlib.pyplot as plt
import numpy as np

from exact_solver import get_hamiltonian, compute_lowest_energies
from vqe_v1 import build_hea

from spsa_optimizer import (
    choose_initial_params,
    run_vqe_spsa_decay_hardware,
    build_and_save_metadata,
    save_summary_and_arrays,
    print_console_summary
)

#---------------------Parameters-------------------------------

# --- Physics ---
N = 4
g = 0.2
N_LAYERS = 1
REAL_ANSATZ = True

# --- SPSA with decay ---
SPSA_MAXITER = 200

A_LR = 0.08
A_SHIFT = 10.0
ALPHA = 0.602
C_PERT = 0.10
GAMMA = 0.101

# --- Shots ---
SHOTS = 8192

# --- Plateau detection ---
PLATEAU_WINDOW = 10
PLATEAU_Z = 1.5
PLATEAU_PATIENCE = 3

# --- Runtime options ---
USE_IDEAL_QPU = False          # Simulates noiselessly while keeping the hardware layout constraint
USE_RESILIENCE_LEVEL_1 = False # set the level of error mitigation but with extra overhead

# --- Initialization ---
INIT_MODE = "resume"       # or "resume" for optimized parameters
INIT_STD_DEV = 0.5

# --- Saving ---
# Force the 'runs' folder to be consistently situated at the project's root directory:
# __file__ is evaluating the absolute path of this python file on ANY running machine.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
RUNS_DIR = PROJECT_ROOT / "runs" / "vqe_fake_marrakesh"
SAVE_QPY = True

# --- Display ---
SHOW_LAYOUT = True
SHOW_CIRCUITS = True

SEED = 42

#main simulation
def main():
    if USE_IDEAL_QPU:
        run_name = f"tfim_vqe_ideal_N{N}_g{g}_L{N_LAYERS}_fake_marrakesh_STAGNANT"
    else:
        run_name = f"tfim_vqe_fake_N{N}_g{g}_L{N_LAYERS}_fake_marrakesh_STAGNANT"
    run_dir = RUNS_DIR / run_name
    run_dir.mkdir(parents=True, exist_ok=True)

    snapshot_npz = run_dir / "latest_snapshot.npz"

    # --- Backend---
    backend = FakeMarrakesh()

    # --- ansatz & hamiltonian ---
    ansatz, _ = build_hea(N, N_LAYERS, real=REAL_ANSATZ)
    n_params = ansatz.num_parameters
    hamiltonian = get_hamiltonian(N, g)

    # --- Transpile ---
    pm = generate_preset_pass_manager(backend=backend, optimization_level=3)
    ansatz_isa = pm.run(ansatz)
    hamiltonian_isa = hamiltonian.apply_layout(layout=ansatz_isa.layout)

    # --- exact values ---
    vals_exact, _ = compute_lowest_energies(N, g)
    E_exact = float(vals_exact[0])

    # --- Save circuits once ---
    if SAVE_QPY:
        with (run_dir / "ansatz_logical.qpy").open("wb") as f:
            qpy.dump(ansatz, f)
        with (run_dir / "ansatz_isa.qpy").open("wb") as f:
            qpy.dump(ansatz_isa, f)

    # --- Display layout & circuits ---
    if ansatz_isa.layout is not None:
        physical_qubits = ansatz_isa.layout.final_index_layout(filter_ancillas=True)
    else:
        physical_qubits = list(range(ansatz_isa.num_qubits))

    if SHOW_LAYOUT:
        try:
            qubit_colors = [
                "#FF4444" if i in physical_qubits else "#DDDDDD"
                for i in range(backend.num_qubits)
            ]
            fig_map = plot_gate_map(backend, qubit_color=qubit_colors, qubit_size=28)
            fig_map.suptitle(f"Qubits used on {backend.name}: {physical_qubits}", fontsize=13)
            fig_map.savefig(run_dir / "qubit_layout.png", dpi=300, bbox_inches="tight")
            plt.close(fig_map)
        except Exception as e:
            print(f"Could not plot gate map: {e}")
            print(f"Physical qubits used: {physical_qubits}")

    if SHOW_CIRCUITS:
        fig_orig = ansatz.draw("mpl")
        if fig_orig:
            fig_orig.savefig(run_dir / "circuit_original.png", dpi=300, bbox_inches="tight")
            plt.close(fig_orig)
        fig_isa = ansatz_isa.draw("mpl", idle_wires=False)
        if fig_isa:
            fig_isa.savefig(run_dir / "circuit_isa.png", dpi=300, bbox_inches="tight")
            plt.close(fig_isa)

    # --- initialization ---
    initial_params = np.array([
        0.71646779, -1.31342932,  0.77267082,  0.96967819, 
       -1.37151579, -1.26823302,  0.48936187, -0.69449997
    ])

    # --- Estimator: same as hardware (job mode for open plan) ---
    if USE_IDEAL_QPU:
        from qiskit_aer.primitives import EstimatorV2 as AerEstimator
        estimator = AerEstimator()
        estimator.options.default_shots = SHOTS
        auth_mode_str = "ideal_qpu_local"
    else:
        estimator = Estimator(mode=backend)
        estimator.options.default_shots = SHOTS
        if USE_RESILIENCE_LEVEL_1:
            estimator.options.resilience_level = 1
        auth_mode_str = "fake_backend_local"

    # --- Metadata saved before first quantum job ---
    build_and_save_metadata(
        run_dir, backend, auth_mode_str,
        N, g, N_LAYERS, REAL_ANSATZ,
        SPSA_MAXITER, A_LR, A_SHIFT, ALPHA, C_PERT, GAMMA,
        SHOTS, USE_RESILIENCE_LEVEL_1,
        PLATEAU_WINDOW, PLATEAU_Z, PLATEAU_PATIENCE,
        INIT_MODE, INIT_STD_DEV, initial_params,
        E_exact, ansatz, ansatz_isa
    )

    print(f"Run directory     : {run_dir}")
    print(f"Backend           : {backend.name}")
    print(f"Ansatz params     : {n_params}")
    print(f"Logical depth     : {ansatz.depth()}")
    print(f"ISA depth         : {ansatz_isa.depth()}")
    print(f"Logical gate count: {dict(ansatz.count_ops())}")
    print(f"ISA gate count    : {dict(ansatz_isa.count_ops())}")
    print(f"Physical qubits   : {physical_qubits}")
    print(f"Exact E0          : {E_exact:.6f}")
    print(f"Shots             : {SHOTS}")
    print(f"Max iterations    : {SPSA_MAXITER}")
    print(f"Ideal QPU         : {USE_IDEAL_QPU}")
    print(f"Resilience        : {USE_RESILIENCE_LEVEL_1}")
    print()

    # --- Run VQE ---
    t0 = time.time()
    result = run_vqe_spsa_decay_hardware(ansatz_isa=ansatz_isa, hamiltonian_isa=hamiltonian_isa,
                                         estimator=estimator, initial_params=initial_params,
                                         run_dir=run_dir,maxiter=SPSA_MAXITER, a_lr=A_LR, a_shift=A_SHIFT, 
                                         alpha=ALPHA,c_pert=C_PERT, gamma=GAMMA,plateau_window=PLATEAU_WINDOW,
                                         plateau_z=PLATEAU_Z, plateau_patience=PLATEAU_PATIENCE,seed=SEED)
    t1 = time.time()

    # --- Save summary ---
    rel_error, iters, early = save_summary_and_arrays(
        run_dir, f"{backend.name} (Ideal QPU)" if USE_IDEAL_QPU else f"{backend.name} (FakeMarrakesh)", result, E_exact, t1 - t0
    )

    # --- Console summary ---
    print_console_summary(
        f"{backend.name} (Ideal QPU)" if USE_IDEAL_QPU else f"{backend.name} (FakeMarrakesh)", E_exact, float(result["energy_opt"]), 
        float(result["energy_opt_se"]), rel_error, iters, early, t1 - t0, run_dir
    )

    # --- Convergence plot ---
    cost_history = result["cost_history"]
    std_history = result["std_history"]
    E_vqe = float(result["energy_opt"])
    
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(1, len(cost_history) + 1)

    ax.plot(x, cost_history, marker="o", markersize=3, linewidth=1.5,
            label=f"E_vqe, err={rel_error*100:.2f}%")
    ax.fill_between(
        x,
        cost_history - std_history,
        cost_history + std_history,
        alpha=0.2,
        label=r"Estimator $\pm 1\sigma$",
    )
    ax.axhline(y=E_exact, linestyle="--", color="red",
               label=f"Exact energy: {E_exact:.6f}")
    ax.axhline(y=E_vqe, linestyle=":", color="green",
               label=f"Reported VQE: {E_vqe:.6f}")

    if early:
        ax.axvline(x=len(cost_history), linestyle=":", alpha=0.6,
                   color="gray", label="Early stop")

    ax.set_xlabel("Iteration")
    ax.set_ylabel(r"Energy $\langle H \rangle$")
    mode_str = "Ideal QPU (Noiseless)" if USE_IDEAL_QPU else "FakeMarrakesh"
    ax.set_title(
        f"VQE on {mode_str} — real HEA, SPSA with decay\n"
        f"(N={N}, g={g}, Layers={N_LAYERS}, Shots={SHOTS})"
    )
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(run_dir / "convergence.png", dpi=300, bbox_inches="tight")
    plt.show()

    print(f"\nConvergence plot saved to {run_dir / 'convergence.png'}")
    print("\nOptimized parameters:")
    print(result["theta_opt"])


if __name__ == "__main__":
    main()
