from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector
import numpy as np
import os
import matplotlib.pyplot as plt
from exact_solver import get_hamiltonian, compute_lowest_energies

#---CONFIGURATION BLOCK---

RUN_DIAGNOSTIC_TEST = False # True = Run Asymptotic Saturation Test instead of standard workflow

NUM_RUNS    = 30      # Number of times to run the script
OPTIMIZER   = "spsa"   # "cobyla" or "spsa"
REAL_ANSATZ = True   # True = Ry-only (2N params), False = Rz-Rx-Rz (6N params)
IDEAL       = True     # True = StatevectorEstimator (no noise), False = noisy FakeMarrakesh QPU simulator

N        = 4  # Number of spins
g        = 0.3  # Transverse field force
N_LAYERS = 2    # Number of HEA repetitions

# --- COBYLA settings ---
COBYLA_MAXITER = 100
COBYLA_RHOBEG  = 0.5
COBYLA_TOL     = 1e-4

# --- SPSA settings ---
SPSA_MAXITER = 222

SHOTS = 4096 # Only used in noisy mode

# --- Quick mode (noisy only) ---
PLOT_LAYOUT_ONLY = False  # True = plot qubit layout and exit (skip VQE)

BASE_SEED = 42

#---ANSATZ CREATION---

def build_hea(N: int, n_layers: int, real: bool = False) -> tuple[QuantumCircuit, ParameterVector]:
    params_per_layer = 2 * N if real else 6 * N
    params = ParameterVector("θ", params_per_layer * n_layers)
    qc = QuantumCircuit(N)

    for layer in range(n_layers):
        offset = params_per_layer * layer

        if real:
            for i in range(N):
                qc.ry(params[offset + i], i)
        else:
            for i in range(N):
                qc.rz(params[offset + i], i)
                qc.rx(params[offset + N + i], i)
                qc.rz(params[offset + 2 * N + i], i)

        qc.barrier()

        for i in range(N):
            qc.cx(i, (i + 1) % N)

        qc.barrier()

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
    cost_history = []

    def cost_func(params):
        pub = (ansatz, [hamiltonian], [params])
        result = estimator.run(pubs=[pub]).result()
        energy = float(result[0].data.evs[0])
        cost_history.append(energy)
        return energy

    from scipy.optimize import minimize
    minimize(cost_func, initial_params, method="cobyla",
             options={"maxiter": cobyla_maxiter,
                      "rhobeg": cobyla_rhobeg,
                      "tol": cobyla_tol})

    E_vqe = np.mean(cost_history[-24:]) if len(cost_history) >= 24 else cost_history[-1]
    return E_vqe, cost_history


#---SPSA VQE---
def run_vqe_qiskit_spsa(ansatz, hamiltonian, estimator, initial_params, maxiter=100):
    cost_history = []

    def cost_func(params):
        pub = (ansatz, [hamiltonian], [params])
        result = estimator.run(pubs=[pub]).result()
        return float(result[0].data.evs[0])
        
    def callback(nfev, x, fx, stepsize, accepted):
        cost_history.append(fx)

    from qiskit_algorithms.optimizers import SPSA
    optimizer = SPSA(maxiter=maxiter, callback=callback)
    
    optimizer.minimize(cost_func, initial_params)
    
    E_vqe = np.mean(cost_history[-24:]) if len(cost_history) >= 24 else (cost_history[-1] if cost_history else 0)
    return E_vqe, cost_history






#------------------------------------------------------------------------------------------------------------------
#---DIAGNOSTIC TEST---
def run_asymptotic_saturation_test(ansatz, hamiltonian, estimator, n_params, initial_seed=42):
    """
    Standalone diagnostic function implementing the Asymptotic Saturation Test
    to determine the optimal number of iterations for SPSA.
    """
    import sys
    M = 100
    SPSA_MAXITER_DIAG = 300
    NOISE_THRESHOLD = 1e-3
    WINDOW = 15
    CONSECUTIVE_REQS = 20
    
    print(f"\n{'='*60}")
    print(f"RUNNING ASYMPTOTIC SATURATION TEST (M={M}, maxiter={SPSA_MAXITER_DIAG})")
    print(f"{'='*60}\n")
    
    all_histories = np.zeros((M, SPSA_MAXITER_DIAG))
    
    for run in range(M):
        np.random.seed(initial_seed + run * 137)
        initial_params = np.random.uniform(-np.pi, np.pi, size=n_params)
        
        # 1. Force optimizer to SPSA
        _, cost_history = run_vqe_qiskit_spsa(ansatz, hamiltonian, estimator, initial_params, maxiter=SPSA_MAXITER_DIAG)
        
        # Populate history
        if len(cost_history) >= SPSA_MAXITER_DIAG:
            all_histories[run, :] = cost_history[:SPSA_MAXITER_DIAG]
        else:
            all_histories[run, :len(cost_history)] = cost_history
            all_histories[run, len(cost_history):] = cost_history[-1]
            
        print(f"  > Diagnostic Run {run+1}/{M} complete.")
        
    print("\n[+] Runs complete. Analyzing convergence...")
    
    # 2. Compute the ensemble mean trajectory
    ensemble_mean = np.mean(all_histories, axis=0)
    
    # 3. Apply a moving tail-average filter (causal window to preserve alignment)
    smoothed_mean = np.zeros_like(ensemble_mean)
    for i in range(len(ensemble_mean)):
        if i < WINDOW:
            smoothed_mean[i] = np.mean(ensemble_mean[:i+1])
        else:
            smoothed_mean[i] = np.mean(ensemble_mean[i-WINDOW+1:i+1])
            
    # 4. Compute absolute discrete derivative
    abs_derivative = np.abs(np.gradient(smoothed_mean))
    
    # 5. Fit exponential decay to the derivative curve: |dE/dk| ≈ a * exp(-b*k) + c
    #    Then solve analytically for k_flat where the fit crosses the noise threshold.
    from scipy.optimize import curve_fit
    
    def exp_decay_deriv(k, a, b, c):
        return a * np.exp(-b * k) + c
    
    x = np.arange(1, SPSA_MAXITER_DIAG + 1)
    
    try:
        p0 = (abs_derivative[0], 0.05, abs_derivative[-1])
        popt, pcov = curve_fit(exp_decay_deriv, x, abs_derivative, p0=p0, maxfev=20000)
        a_fit, b_fit, c_fit = popt
        fit_curve = exp_decay_deriv(x, *popt)
        
        # R² of the fit
        ss_res = np.sum((abs_derivative - fit_curve) ** 2)
        ss_tot = np.sum((abs_derivative - np.mean(abs_derivative)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        
        # Analytical k_flat: solve a*exp(-b*k) + c = epsilon  =>  k = ln(a / (epsilon - c)) / b
        if NOISE_THRESHOLD > c_fit and a_fit > 0 and b_fit > 0:
            k_flat = int(np.log(a_fit / (NOISE_THRESHOLD - c_fit)) / b_fit)
        else:
            # Fallback: if c >= threshold (floor is already above threshold), use where fit is 2x floor
            k_flat = int(np.log(a_fit / c_fit) / b_fit) if b_fit > 0 and a_fit > 0 and c_fit > 0 else SPSA_MAXITER_DIAG
            
        k_flat = max(1, min(k_flat, SPSA_MAXITER_DIAG))
        k_optimal = int(k_flat * 1.2)
        
        print(f"  -> Exponential Fit: |dE/dk| ≈ {a_fit:.4f} * exp(-{b_fit:.5f} * k) + {c_fit:.6f}")
        print(f"  -> Fit R²                            : {r_squared:.5f}")
        print(f"  -> Flattening Detection (k_flat)     : {k_flat} iterations")
        print(f"  -> Recommended SPSA_MAXITER (k_opt)  : {k_optimal} iterations")
        fit_success = True
        
    except Exception as e:
        print(f"  WARNING: Exponential fit failed ({e}). Falling back to threshold crossing.")
        k_flat = SPSA_MAXITER_DIAG
        k_optimal = SPSA_MAXITER_DIAG
        fit_success = False
        
    # 6. Visualization
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 8))
    
    # Top Subplot
    ax1.plot(x, ensemble_mean, color='gray', alpha=0.5, label='Raw ensemble mean')
    ax1.plot(x, smoothed_mean, color='#0066CC', linewidth=2, label=f'Smoothed trajectory (window={WINDOW})')
    if k_flat < SPSA_MAXITER_DIAG:
        ax1.axvline(x=k_flat, color='red', linestyle='--', linewidth=2, label=f'k_flat = {k_flat}')
    ax1.set_ylabel(r"Energy $\langle H \rangle$")
    ax1.set_title(f"Asymptotic saturation test ({M} runs)")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Bottom Subplot
    ax2.plot(x, abs_derivative, color='purple', linewidth=1.5, alpha=0.6, label=r'$|d\langle E \rangle / dk|$')
    if fit_success:
        ax2.plot(x, fit_curve, color='orange', linewidth=2, 
                 label=f'Exp fit ($R^2$={r_squared:.4f})')
    ax2.axhline(y=NOISE_THRESHOLD, color='black', linestyle=':', label=f'Noise threshold ({NOISE_THRESHOLD})')
    if k_flat < SPSA_MAXITER_DIAG:
        ax2.axvline(x=k_flat, color='red', linestyle='--', linewidth=2, label=f'k_flat = {k_flat}')
    ax2.set_ylabel("Absolute derivative")
    ax2.set_xlabel("Optimizer iterations")
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig("asymptotic_saturation_test.png", dpi=300)
    print(f"  -> Diagnostic plot saved to asymptotic_saturation_test.png")
    plt.show()
    sys.exit(0)

#------------------------------------------------------------------------------------------

if __name__ == "__main__":

    tag = f"{'real' if REAL_ANSATZ else 'general'}_{'ideal' if IDEAL else 'noisy'}_{OPTIMIZER}_avg{NUM_RUNS}"
    FIG_DIR = f"figures/vqe_figs/{tag}"
    os.makedirs(FIG_DIR, exist_ok=True)

    ansatz, params = build_hea(N, N_LAYERS, real=REAL_ANSATZ)
    n_params = ansatz.num_parameters
    Hamiltonian = get_hamiltonian(N, g)

    if IDEAL:
        from qiskit.primitives import StatevectorEstimator
        print(f"MODE: IDEAL")
        estimator = StatevectorEstimator()
        circuit_to_run = ansatz
        hamiltonian_to_run = Hamiltonian
        shots_label = "ideal"
    else:
        from qiskit_ibm_runtime import EstimatorV2 as Estimator, Batch
        from qiskit_ibm_runtime.fake_provider import FakeMarrakesh
        from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

        backend = FakeMarrakesh()
        pm = generate_preset_pass_manager(target=backend.target, optimization_level=3)
        ansatz_isa = pm.run(ansatz)
        Hamiltonian_isa = Hamiltonian.apply_layout(layout=ansatz_isa.layout)

        circuit_to_run = ansatz_isa
        hamiltonian_to_run = Hamiltonian_isa
        shots_label = f"{SHOTS}shots"

    # --- Run Diagnostic Test Intercept ---
    if RUN_DIAGNOSTIC_TEST:
        # Pass the constructed estimator/circuit matching the IDEAL/Noisy mode to ensure consistency
        if not IDEAL:
            # Recreate an unbatched estimator just for the diagnostic
            estimator_diag = Estimator(mode=backend)
            estimator_diag.options.default_shots = SHOTS
        else:
            estimator_diag = estimator

        run_asymptotic_saturation_test(circuit_to_run, hamiltonian_to_run, estimator_diag, n_params, initial_seed=BASE_SEED)

    vals_exact, _ = compute_lowest_energies(N, g)
    E_exact = float(vals_exact[0])

    all_cost_histories = []
    all_E_vqes = []

    print(f"\n{'='*50}")
    print(f"Starting {NUM_RUNS} runs of VQE...")
    print(f"Exact energy: {E_exact:.6f}")
    print(f"{'='*50}\n")

    if not IDEAL and OPTIMIZER == "cobyla":
        # Keep batch open for COBYLA noisy runs
        batch = Batch(backend=backend)
        estimator = Estimator(mode=batch)
        estimator.options.default_shots = SHOTS

    if not IDEAL and OPTIMIZER == "spsa":
        estimator = Estimator(mode=backend)
        estimator.options.default_shots = SHOTS

    for run_idx in range(NUM_RUNS):
        # New seed for each run
        np.random.seed(BASE_SEED + run_idx)
        initial_params = np.random.uniform(-np.pi, np.pi, size=n_params)

        print(f"--- Run {run_idx + 1}/{NUM_RUNS} ---")

        if OPTIMIZER == "cobyla":
            E_vqe, cost_history = run_vqe_cobyla(
                circuit_to_run, hamiltonian_to_run, estimator, initial_params,
                cobyla_maxiter=COBYLA_MAXITER, cobyla_rhobeg=COBYLA_RHOBEG, cobyla_tol=COBYLA_TOL
            )
        else:
            E_vqe, cost_history = run_vqe_qiskit_spsa(
                ansatz=circuit_to_run, hamiltonian=hamiltonian_to_run, estimator=estimator,
                initial_params=initial_params, maxiter=SPSA_MAXITER
            )
        
        all_cost_histories.append(cost_history)
        all_E_vqes.append(E_vqe)
        print(f"    E_vqe = {E_vqe:.6f} (Error: {abs(E_vqe - E_exact)/abs(E_exact)*100:.2f}%)")

    if not IDEAL and OPTIMIZER == "cobyla":
        batch.close()

    #---PROCESS AND PLOT RESULTS---

    # Pad cost histories to max length if they stopped at different times
    max_len = max(len(h) for h in all_cost_histories)
    padded_histories = []
    for h in all_cost_histories:
        pad_len = max_len - len(h)
        padded_histories.append(np.pad(h, (0, pad_len), 'edge'))
    
    all_histories = np.array(padded_histories)
    
    mean_history = np.mean(all_histories, axis=0)
    median_history = np.median(all_histories, axis=0)
    

    # Polyak-Ruppert Averaging: Mean of the last 37 iterations for each run
    # Window derived from saturation test: [k_flat=185, K_max=222] => 37 iterations
    tail_averages = np.mean(all_histories[:, -37:], axis=1)
    best_run_idx = np.argmin(tail_averages)
    E_vqe = tail_averages[best_run_idx]  # Minimal tail average
    
    mean_E_vqe = np.mean(tail_averages)
    mean_rel_error = abs(mean_E_vqe - E_exact) / abs(E_exact)
    best_rel_error = abs(E_vqe - E_exact) / abs(E_exact)

    print(f"\n{'─'*56}")
    print(f"  Runs finalized : {NUM_RUNS}")
    print(f"  Exact energy   : {E_exact:.6f}")
    print(f"  Mean VQE energy: {mean_E_vqe:.6f}")
    print(f"  Best VQE (avg) : {E_vqe:.6f} (Error: {best_rel_error*100:.3f}%)")
    print(f"{'─'*56}")


    # -------------------------------------Plot-------------------------------------


    ansatz_label = "real HEA" if REAL_ANSATZ else "general HEA"
    sim_label = "ideal" if IDEAL else "noisy"
    optim_label = OPTIMIZER.upper()

    x = np.arange(1, max_len + 1)

    # --- Main Plot ---
    plt.figure(figsize=(10, 5))

    # Identify the best run index already calculated via tail average
    best_run_history = all_histories[best_run_idx]

    # Plot all runs in gray with high transparency
    for h in all_histories:
        plt.plot(x, h, color='gray', alpha=0.15, linewidth=1)

    # Plot the median in dashed black
    plt.plot(x, median_history, color='black', linestyle='--', linewidth=2, 
             label="Median (all runs)")

    # Best run trajectory
    plt.plot(x, best_run_history, color='#00CC66', linewidth=2.5, zorder=10,
             label=f"Best run (err={best_rel_error*100:.2f}%)")

    # Exact energy
    plt.axhline(y=E_exact, color='red', linestyle='--', linewidth=2, zorder=5,
                label=f"Exact energy: {E_exact:.6f}")

    plt.xlabel("Optimizer iterations")
    plt.ylabel(r"Energy $\langle H \rangle$")
    plt.title(f"Mean energy over {NUM_RUNS} runs, {optim_label} on {sim_label} QPU\n"
              f"(N={N}, g={g}, Layers={N_LAYERS})")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plot_path = f"{FIG_DIR}/convergence_avg_{NUM_RUNS}_N{N}_g{g}_{shots_label}.png"
    plt.savefig(plot_path, dpi=300)
    print(f"\nConvergence plot saved to {plot_path}")

    # --- Second Plot: Exponential Fit ---
    from scipy.optimize import curve_fit
    def exp_decay(t, a, b, c):
        return a * np.exp(-b * t) + c
        
    try:
        p0 = (mean_history[0] - mean_history[-1], 0.1, mean_history[-1])
        popt, pcov = curve_fit(exp_decay, x, mean_history, p0=p0, maxfev=10000)
        fit_y = exp_decay(x, *popt)
        
        # Calculate R^2
        ss_res = np.sum((mean_history - fit_y) ** 2)
        ss_tot = np.sum((mean_history - np.mean(mean_history)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        
        plt.figure(figsize=(10, 5))
    
        plt.plot(x, mean_history, 'o', color='#0066CC', alpha=0.5, markersize=5, label="Mean energy")
        
        # Plot Fit
        fit_label = f"Exp fit: ${popt[0]:.3f}e^{{{-popt[1]:.4f}t}} {popt[2]:+.3f}$\n$R^2 = {r_squared:.5f}$"
        plt.plot(x, fit_y, color='purple', linestyle='-', linewidth=2.5, label=fit_label)
        
        plt.xlabel("Optimizer iterations")
        plt.ylabel(r"Mean energy $\langle H \rangle$")
        plt.title(f"Exponential fit of mean VQE convergence ({NUM_RUNS} runs)\n"
                  f"(N={N}, g={g}, Layers={N_LAYERS})")
        
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)
        fit_path = f"{FIG_DIR}/convergence_avg_{NUM_RUNS}_N{N}_g{g}_{shots_label}_FIT.png"
        plt.savefig(fit_path, dpi=300)
        print(f"Fit plot saved to {fit_path}")
    except Exception as e:
        print(f"Could not calculate exponential fit: {e}")

    plt.show()
