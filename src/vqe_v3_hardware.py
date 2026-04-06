"""
Production VQE on a REAL IBM QPU (Open Plan friendly, no session).

Purpose
-------
- Run TFIM VQE on a real IBM backend in job mode (no Session).
- Save all useful information incrementally, so a partially completed run is still usable.
- Do NOT plot during execution.
- Use manual SPSA with decaying gain sequences and noise-aware plateau detection.

Authentication
--------------
Option 1 (recommended for this script):
    export IQP_TOKEN="your_api_key"
    export IBM_INSTANCE_CRN="your_instance_crn"   # optional, only if needed
    python vqe_hardware_production.py

Option 2 (trusted environment):
    from qiskit_ibm_runtime import QiskitRuntimeService
    QiskitRuntimeService.save_account(token="...", instance="...", overwrite=True)
    # then the script can use QiskitRuntimeService() without env vars
"""

from __future__ import annotations

import os
import time
import datetime as dt
from pathlib import Path

from qiskit import qpy
from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2 as Estimator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

from exact_solver import get_hamiltonian, compute_lowest_energies
from vqe_v1 import build_hea

from spsa_optimizer import (
    choose_initial_params,
    run_vqe_spsa_decay_hardware,
    build_and_save_metadata,
    save_summary_and_arrays,
    print_console_summary
)

# ═══════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════

BACKEND_NAME = "ibm_marrakesh"

# --- Physics ---
N = 4
g = 0.2
N_LAYERS = 1
REAL_ANSATZ = True

# --- SPSA with decay ---
SPSA_MAXITER = 30

# ak = A_LR / (k + 1 + A_SHIFT)^ALPHA
A_LR = 0.08
A_SHIFT = 10.0
ALPHA = 0.602

# ck = C_PERT / (k + 1)^GAMMA
C_PERT = 0.10
GAMMA = 0.101

# --- Shots ---
SHOTS = 4096

# --- Plateau detection ---
PLATEAU_WINDOW = 10
PLATEAU_Z = 1.5
PLATEAU_PATIENCE = 3

# --- Runtime options ---
USE_RESILIENCE_LEVEL_1 = False   # set True if you want simple mitigation, at extra overhead

# --- Initialization ---
INIT_MODE = "random"          # "random" | "resume"
INIT_STD_DEV = 0.05

# --- Saving ---
# Force the 'runs' folder to be consistently situated at the project's root directory:
# __file__ is evaluating the absolute path of this python file on ANY running machine.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
RUNS_DIR = PROJECT_ROOT / "runs" / "vqe_hardware"
SAVE_QPY = True

SEED = 42

# ═══════════════════════════════════════════════════════════════════

RUNS_DIR.mkdir(parents=True, exist_ok=True)

def main():
    run_name = f"tfim_vqe_hw_N{N}_g{g}_L{N_LAYERS}_{BACKEND_NAME}"
    run_dir = RUNS_DIR / run_name
    run_dir.mkdir(parents=True, exist_ok=True)

    snapshot_npz = run_dir / "latest_snapshot.npz"

    # --- Authentication / service ---
    token = os.environ.get("IQP_TOKEN", None)
    instance = os.environ.get("IBM_INSTANCE_CRN", None)

    if token:
        service = QiskitRuntimeService(token=token, instance=instance)
        auth_mode = "env_token"
    else:
        service = QiskitRuntimeService(instance=instance)
        auth_mode = "saved_credentials"

    backend = service.backend(BACKEND_NAME)

    # --- Build ansatz / Hamiltonian ---
    ansatz, _ = build_hea(N, N_LAYERS, real=REAL_ANSATZ)
    n_params = ansatz.num_parameters
    hamiltonian = get_hamiltonian(N, g)

    # --- Transpile ---
    pm = generate_preset_pass_manager(backend=backend, optimization_level=3)
    ansatz_isa = pm.run(ansatz)
    hamiltonian_isa = hamiltonian.apply_layout(layout=ansatz_isa.layout)

    # --- Exact reference ---
    vals_exact, _ = compute_lowest_energies(N, g)
    E_exact = float(vals_exact[0])

    # --- Save circuits once ---
    if SAVE_QPY:
        with (run_dir / "ansatz_logical.qpy").open("wb") as f:
            qpy.dump(ansatz, f)
        with (run_dir / "ansatz_isa.qpy").open("wb") as f:
            qpy.dump(ansatz_isa, f)

    # --- Choose init ---
    initial_params = choose_initial_params(
        n_params=n_params,
        snapshot_npz=snapshot_npz,
        init_mode=INIT_MODE,
        init_std_dev=INIT_STD_DEV,
        seed=SEED,
    )

    # --- Estimator ---
    estimator = Estimator(mode=backend)
    estimator.options.default_shots = SHOTS
    if USE_RESILIENCE_LEVEL_1:
        estimator.options.resilience_level = 1

    # --- Metadata saved before first quantum job ---
    build_and_save_metadata(
        run_dir, backend, auth_mode,
        N, g, N_LAYERS, REAL_ANSATZ,
        SPSA_MAXITER, A_LR, A_SHIFT, ALPHA, C_PERT, GAMMA,
        SHOTS, USE_RESILIENCE_LEVEL_1,
        PLATEAU_WINDOW, PLATEAU_Z, PLATEAU_PATIENCE,
        INIT_MODE, INIT_STD_DEV, initial_params,
        E_exact, ansatz, ansatz_isa
    )

    print(f"Run directory   : {run_dir}")
    print(f"Backend         : {backend.name}")
    print(f"Ansatz params   : {n_params}")
    print(f"ISA depth       : {ansatz_isa.depth()}")
    print(f"ISA gate count  : {dict(ansatz_isa.count_ops())}")
    print(f"Exact E0        : {E_exact:.6f}")
    print(f"Shots           : {SHOTS}")
    print(f"Max iterations  : {SPSA_MAXITER}")
    print(f"Resilience L1   : {USE_RESILIENCE_LEVEL_1}")
    print()

    # --- Run ---
    t0 = time.time()
    result = run_vqe_spsa_decay_hardware(
        ansatz_isa=ansatz_isa,
        hamiltonian_isa=hamiltonian_isa,
        estimator=estimator,
        initial_params=initial_params,
        run_dir=run_dir,
        maxiter=SPSA_MAXITER,
        a_lr=A_LR,
        a_shift=A_SHIFT,
        alpha=ALPHA,
        c_pert=C_PERT,
        gamma=GAMMA,
        plateau_window=PLATEAU_WINDOW,
        plateau_z=PLATEAU_Z,
        plateau_patience=PLATEAU_PATIENCE,
        seed=SEED,
    )
    t1 = time.time()

    rel_error, iters, early = save_summary_and_arrays(
        run_dir, backend.name, result, E_exact, t1 - t0
    )

    print_console_summary(
        backend.name, E_exact, float(result["energy_opt"]), 
        float(result["energy_opt_se"]), rel_error, iters, early, t1 - t0, run_dir
    )


if __name__ == "__main__":
    main()