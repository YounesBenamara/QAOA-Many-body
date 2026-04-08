
# Variational Ising QPT

This project combines classical simulation and quantum hardware execution to implement a variational quantum-classical framework for detecting and characterizing the quantum phase transition in the one-dimensional transverse-field Ising model (TFIM). The repository studies how the ground state evolves as the transverse-field parameter $g$ is varied across the critical point $g = 1$, where the system transitions from a ferromagnetic ordered phase to a paramagnetic disordered phase.

## Current Status

This project is currently under development. At present:

- The implementation of the Quantum Approximate Optimization Algorithm (QAOA) is not yet complete.
- A formal bibliography and explicit source citations have not yet been integrated into the documentation.
- The repository layout and file organization may still change as the project evolves.
- The `figures/` directory is being reorganized for clarity and consistency.
- A detailed report presenting the full results and physical analysis is being prepared in LaTeX.

## Project Overview

The repository explores the ability of noisy intermediate-scale quantum (NISQ) algorithms to capture many-body physics signatures using three complementary approaches:

1. **Exact diagonalization**: a reference baseline using SciPy to compute the exact ground-state energy and physical observables.
2. **Variational Quantum Eigensolver (VQE)**: a variational ground-state solver based on a hardware-efficient ansatz optimized with SPSA and COBYLA.
3. **Quantum Approximate Optimization Algorithm (QAOA)**: a structured variational ansatz built directly from the Hamiltonian terms, intended to probe phase-sensitive observables.

## Physical Model

The system considered is a ring of $N$ spin- $\tfrac{1}{2}$ particles with periodic boundary conditions, governed by the Hamiltonian

$$
H(g) = - \sum_{i=0}^{N-1} Z_i Z_{i+1} - g \sum_{i=0}^{N-1} X_i
$$

where $Z_i$ and $X_i$ are Pauli operators acting on site $i$, and periodicity implies that $Z_N \equiv Z_0$.

The model exhibits two main regimes:

- **Ordered phase ($g < 1$)**: the nearest-neighbor interaction term dominates, favoring ferromagnetic ordering.
- **Disordered phase ($g > 1$)**: the transverse field dominates, and the ground state tends toward a product state polarized along the $x$-direction.

At the critical point $g = 1$, the model undergoes a quantum phase transition in the thermodynamic limit. For finite systems, such as $N = 4$, one observes a finite-size crossover together with precursors of critical behavior rather than a true non-analytic phase transition.

## Repository Structure

```bash
Variational_Ising_QPT/
├── runs/
│   ├── vqe_fake_marrakesh/
│   └── vqe_hardware/
├── src/
│   ├── figures/
│   │   ├── exact_solver_figs/
│   │   └── vqe_figs/
│   ├── exact_solver.py
│   ├── g_comparison.py
│   ├── layers_comparison.py
│   ├── saturation.md
│   ├── spsa_optimizer.py
│   ├── vqe_v1.py
│   ├── vqe_v2.py
│   └── vqe_v3_hardware.py
└── README.md
```

## Getting Started

### Prerequisites

A Python environment with the following packages is required:

- `qiskit`
- `qiskit-ibm-runtime`
- `qiskit-algorithms`
- `numpy`
- `scipy`
- `matplotlib`

It is recommended to install dependencies in a dedicated virtual environment.

### Running the Diagnostic

To estimate the iteration budget needed for a given simulator or hardware backend, the asymptotic saturation diagnostic implemented in the VQE scripts can be used:

1. Set the diagnostic flag to `True` in the configuration block of the relevant source file.
2. Run the script.
3. Analyze the generated convergence curve and exponential fit to estimate the onset of saturation.

## Key Results Summary

- **Optimal depth**: for a four-qubit system, a two-layer hardware-efficient ansatz appears to provide a good compromise between expressivity and trainability.
- **Convergence**: in ideal simulations, SPSA with decaying step sizes can reach energies within about 1% of the exact diagonalization reference.
- **Phase signatures**: the framework tracks physically relevant observables, including the suppression of the magnetic structure factor and the peak of bipartite entanglement entropy near the critical region.

## Notes

This repository is intended both as a computational study of the transverse-field Ising model and as a practical investigation of how well current variational quantum algorithms can recover many-body signatures on simulators and real quantum hardware.