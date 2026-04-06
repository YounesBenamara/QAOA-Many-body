
from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector
import numpy as np
import os
import matplotlib.pyplot as plt
from exact_solver import get_hamiltonian, compute_lowest_energies


#---CONFIGURATION BLOCK---

OPTIMIZER   = "spsa"   # "cobyla" or "spsa"
REAL_ANSATZ = True     # True = Ry-only (2N params), False = Rz-Rx-Rz (6N params)
IDEAL       = True  # True = StatevectorEstimator (no noise), False = noisy FakeMarrakesh QPU simulator

N        = 4  # Number of spins
g        = 0.3  # Transverse field force
N_LAYERS = 4    # Number of HEA repetitions

# --- COBYLA settings ---
COBYLA_MAXITER = 100
COBYLA_RHOBEG  = 0.5
COBYLA_TOL     = 1e-4

# --- SPSA settings ---
SPSA_MAXITER = 150

SHOTS = 4096 # Only used in noisy mode

# --- Quick mode (noisy only) ---
PLOT_LAYOUT_ONLY = False  # True = plot qubit layout and exit (skip VQE)

SEED = 42
np.random.seed(SEED)

#---ANSATZ CREATION---


def build_hea(N: int, n_layers: int, real: bool = False) -> tuple[QuantumCircuit, ParameterVector]:
    '''
    Hardware-Efficient Ansatz (HEA) for the TFIM.

    If real=False (general):
        Each layer: Rz-Rx-Rz rotations → CNOT → Rz-Rx-Rz rotations
        → 6N parameters per layer.
    If real=True (real-valued):
        Each layer: Ry rotations → CNOT → Ry rotations
        2N parameters per layer.

    Args:
        N        : number of qubits
        n_layers : number of HEA repetitions
        real     : if True, use Ry-only (real-valued ansatz)

    Returns:
        (circuit, params): parameterized QuantumCircuit and ParameterVector
    '''
    params_per_layer = 2 * N if real else 6 * N
    params = ParameterVector("θ", params_per_layer * n_layers)
    qc = QuantumCircuit(N)

    for layer in range(n_layers):
        offset = params_per_layer * layer

        # Pre-entanglement rotations
        if real:
            for i in range(N):
                qc.ry(params[offset + i], i)
        else:
            for i in range(N):
                qc.rz(params[offset + i], i)
                qc.rx(params[offset + N + i], i)
                qc.rz(params[offset + 2 * N + i], i)

        qc.barrier()

        # Entangling layer (periodic boundary conditions)
        for i in range(N):
            qc.cx(i, (i + 1) % N)

        qc.barrier()

        # Post-entanglement rotations
        if real:
            for i in range(N):
                qc.ry(params[offset + N + i], i)
        else:
            for i in range(N):
                qc.rz(params[offset + 3 * N + i], i)
                qc.rx(params[offset + 4 * N + i], i)
                qc.rz(params[offset + 5 * N + i], i)

        qc.barrier()

    return qc, params


#---COBYLA VQE---
def run_vqe_cobyla(ansatz, hamiltonian, estimator, initial_params,
                   cobyla_maxiter=50, cobyla_rhobeg=0.5, cobyla_tol=1e-4):
    '''
    Run VQE with COBYLA optimizer.

    Returns:
        (E_vqe, cost_history)
    '''
    cost_history = []

    def cost_func(params):
        pub = (ansatz, [hamiltonian], [params])
        result = estimator.run(pubs=[pub]).result()
        energy = result[0].data.evs[0]
        cost_history.append(energy)
        if not IDEAL:
            print(f"Iteration {len(cost_history)} - Energy: {energy}")
        return energy

    from scipy.optimize import minimize
    result = minimize(cost_func, initial_params, method="cobyla",
                      options={"maxiter": cobyla_maxiter,
                               "rhobeg": cobyla_rhobeg,
                               "tol": cobyla_tol})

    # Average the last 24 entries (Polyak-Ruppert window from saturation test)
    E_vqe = np.mean(cost_history[-24:])

    return E_vqe, cost_history


#---SPSA VQE---
def run_vqe_qiskit_spsa(ansatz, hamiltonian, estimator, initial_params, maxiter=100):
    '''
    Run VQE using Qiskit's built-in SPSA optimizer algorithm.
    '''
    cost_history = []

    def cost_func(params):
        pub = (ansatz, [hamiltonian], [params])
        result = estimator.run(pubs=[pub]).result()
        return float(result[0].data.evs[0])
        
    def callback(nfev, x, fx, stepsize, accepted):
        cost_history.append(fx)
        if not IDEAL:
            print(f"Iteration {len(cost_history)} - Energy: {fx}")

    from qiskit_algorithms.optimizers import SPSA
    optimizer = SPSA(maxiter=maxiter, callback=callback)
    
    result = optimizer.minimize(cost_func, initial_params)
    
    # Polyak-Ruppert averaging over last 24 iterations (from saturation test)
    E_vqe = np.mean(cost_history[-24:])

    return E_vqe, cost_history

#------------------------------------------------------------------------------------------

if __name__ == "__main__":

    # --- Output directories ---
    tag = f"{'real' if REAL_ANSATZ else 'general'}_{'ideal' if IDEAL else 'noisy'}_{OPTIMIZER}"
    FIG_DIR = f"figures/vqe_figs/{tag}"
    os.makedirs(FIG_DIR, exist_ok=True)

    if not IDEAL:
        CHECKPOINT = f"{FIG_DIR}/params_N{N}_g{g}.npy"

    #---ANSATZ & BACKEND SETUP---

    # Build ansatz
    ansatz, params = build_hea(N, N_LAYERS, real=REAL_ANSATZ)
    n_params = ansatz.num_parameters

    # Hamiltonian
    Hamiltonian = get_hamiltonian(N, g)

    if IDEAL:
        # --- Ideal simulation: no noise, no transpilation ---
        from qiskit.primitives import StatevectorEstimator

        print(f"MODE: IDEAL")
        print(f"Ansatz depth      : {ansatz.depth()}")
        print(f"Ansatz gate count : {dict(ansatz.count_ops())}")
        print(f"Number of params  : {n_params}")

        estimator = StatevectorEstimator()
        circuit_to_run = ansatz
        hamiltonian_to_run = Hamiltonian
        shots_label = "ideal"

    else:
        # --- Noisy simulation on fake backend (FakeMarrakesh) ---

        from qiskit_ibm_runtime import EstimatorV2 as Estimator, Batch
        from qiskit_ibm_runtime.fake_provider import FakeMarrakesh
        from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

        backend = FakeMarrakesh()
        pm = generate_preset_pass_manager(target=backend.target, optimization_level=3)

        ansatz_isa = pm.run(ansatz)

        print(f"Mode: Fakebackend simulation")
        print(f"Original depth      : {ansatz.depth()}")
        print(f"ISA depth           : {ansatz_isa.depth()}")
        print(f"Original gate count : {dict(ansatz.count_ops())}")
        print(f"ISA gate count      : {dict(ansatz_isa.count_ops())}")

        Hamiltonian_isa = Hamiltonian.apply_layout(layout=ansatz_isa.layout)

        # --- Plot physical qubit layout on QPU ---
        try:
            from qiskit.visualization import plot_gate_map

            physical_qubits = ansatz_isa.layout.final_index_layout(filter_ancillas=True)
            print(f"Physical qubits : {physical_qubits}")

            qubit_colors = ['blue' if i in physical_qubits else 'gray'
                            for i in range(backend.num_qubits)]

            fig_map = plot_gate_map(backend, qubit_color=qubit_colors,
                                    qubit_size=28)

            fig_map.suptitle(f"Qubits used on {backend.name}: {physical_qubits}", fontsize=13)
            fig_map.savefig(f"{FIG_DIR}/qubit_layout_N{N}_g{g}.png", dpi=300, bbox_inches='tight')
            print(f"Qubit layout saved to {FIG_DIR}/qubit_layout_N{N}_g{g}.png")
        except Exception as e:
            print(f"Could not plot gate map: {e}")

        # --- Plot circuits ---
        #Ideal circuit
        fig_circ = ansatz.draw('mpl')
        fig_circ.savefig(f"{FIG_DIR}/circuit_original_N{N}_g{g}.png", dpi=300, bbox_inches='tight')
        print(f"Original circuit saved to {FIG_DIR}/circuit_original_N{N}_g{g}.png")

        #ISA circuit
        fig_isa = ansatz_isa.draw('mpl', idle_wires=False)
        fig_isa.savefig(f"{FIG_DIR}/circuit_isa_N{N}_g{g}.png", dpi=300, bbox_inches='tight')
        print(f"ISA circuit saved to {FIG_DIR}/circuit_isa_N{N}_g{g}.png")

        plt.show()

        if PLOT_LAYOUT_ONLY:
            raise SystemExit(0)

        circuit_to_run = ansatz_isa
        hamiltonian_to_run = Hamiltonian_isa
        shots_label = f"{SHOTS}shots"

    #---Exact energy---
    vals_exact, _ = compute_lowest_energies(N, g)
    E_exact = float(vals_exact[0])

    #---Initial parameters---

    if not IDEAL and os.path.exists(CHECKPOINT):
        #starting from last computed parameters for noisy simulation
        initial_params = np.load(CHECKPOINT)
        print(f"Loaded checkpoint params from {CHECKPOINT}")
    else:
        # Generic non-optimized random start Uniform[-pi, pi]
        initial_params = np.random.uniform(-np.pi, np.pi, size=n_params)
        print("Starting from totally random parameters uniformly in [-pi, pi]")


    #---OPTIMIZATION---

    if OPTIMIZER == "cobyla":
        # --- COBYLA ---
        if not IDEAL:
            batch = Batch(backend=backend)
            estimator = Estimator(mode=batch)
            estimator.options.default_shots = SHOTS

        E_vqe, cost_history = run_vqe_cobyla(
            circuit_to_run, hamiltonian_to_run, estimator, initial_params,
            cobyla_maxiter=COBYLA_MAXITER, cobyla_rhobeg=COBYLA_RHOBEG, cobyla_tol=COBYLA_TOL
        )
        std_history = None

        if not IDEAL:
            batch.close()

    elif OPTIMIZER == "spsa":
        # --- Generic SPSA ---
        if not IDEAL:
            estimator = Estimator(mode=backend)
            estimator.options.default_shots = SHOTS

        E_vqe, cost_history = run_vqe_qiskit_spsa(
            ansatz=circuit_to_run,
            hamiltonian=hamiltonian_to_run,
            estimator=estimator,
            initial_params=initial_params,
            maxiter=SPSA_MAXITER,
        )
        std_history = None


    #---RESULTS---

    rel_error = abs(E_vqe - E_exact) / abs(E_exact)

    print(f"\n{'─'*56}")
    print(f"  Config         : {tag}")
    print(f"  Exact energy   : {E_exact:.6f}")
    if OPTIMIZER == "spsa" and not IDEAL:
        E_vqe_se = result["energy_opt_se"]
        print(f"  VQE energy     : {E_vqe:.6f} ± {E_vqe_se:.6f}")
    else:
        print(f"  VQE energy     : {E_vqe:.6f}")
    print(f"  Relative error : {rel_error*100:.3f} %")
    print(f"  Iterations     : {len(cost_history)}")
    if std_history is not None:
        print(f"  Max std        : {np.max(std_history):.6f}")
    print(f"{'─'*56}")


    #---CONVERGENCE PLOT---

    ansatz_label = "real HEA" if REAL_ANSATZ else "general HEA"
    sim_label = "ideal" if IDEAL else "noisy"
    optim_label = OPTIMIZER.upper()

    plt.figure(figsize=(10, 5))
    x = np.arange(1, len(cost_history) + 1)

    plt.plot(x, cost_history, marker="o", markersize=3, linewidth=1.5,
             label=f"E_vqe, err={rel_error*100:.2f}%", color='#0066CC')

    if std_history is not None and not IDEAL:
        plt.fill_between(
            x,
            cost_history - std_history,
            cost_history + std_history,
            alpha=0.2,
            label=r"Estimator $\pm 1\sigma$",
        )

    plt.axhline(y=E_exact, color='#FF3300', linestyle='--',
                label=f"Exact Energy: {E_exact:.6f}")
    plt.axhline(y=E_vqe, color='#00AA44', linestyle=':',
                label=f"Reported VQE: {E_vqe:.6f}")

    plt.xlabel("Iteration")
    plt.ylabel(r"Energy $\langle H \rangle$")
    plt.title(f"VQE with {ansatz_label}, {optim_label} on {sim_label}"
              f" (N={N}, g={g}, Layers={N_LAYERS})")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{FIG_DIR}/convergence_N{N}_g{g}_{shots_label}.png", dpi=300)
    plt.show()
