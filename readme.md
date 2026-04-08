# Variational Ising QPT

This project utilizes both classical simulation and real hardware execution to implement a variational quantum-classical framework designed to detect and characterize the quantum phase transition in the one-dimensional transverse-field Ising model. By utilizing variational algorithms, the repository tracks the evolution of the ground state across the critical point g = 1, transitioning from a ferromagnetic ordered phase to a paramagnetic disordered phase.

## Current Status

This project is currently under development and the following conditions apply:
* The implementation of the Quantum Approximate Optimization Algorithm is currently lacking.
* A comprehensive list of formal sources and citations is yet to be incorporated into the final documentation.
* The repository layout and file organization are subject to change as the project evolves.
* The figures directory is currently being reorganized and cleaned for better clarity.
* A detailed report documenting the full findings and physical analysis is currently in preparation using LaTeX.

## Project Overview

[span_0](start_span)The repository explores the capabilities of algorithms from the noisy intermediate-scale quantum era to capture many-body physics signatures using three complementary approaches[span_0](end_span):

1. Exact Diagonalization: A reference baseline utilizing SciPy for benchmarking energy and physical observables.
2. Variational Quantum Eigensolver: Ground state approximation using a hardware-efficient ansatz with optimized SPSA and COBYLA loops.
3. QAOA: A structured variational strategy tied directly to the Hamiltonian terms intended to probe phase-sensitive observables.

## Physical Model

The simulation considers a ring of N spin-1/2 particles with periodic boundary conditions governed by the following Hamiltonian:

$$H(g) = - \sum_{i=0}^{N-1} Z_i Z_{i+1} - g \sum_{i=0}^{N-1} X_i$$

The model exhibits two primary regimes:
* Ordered Phase (g < 1): This regime is dominated by nearest-neighbor interactions where the ground state exhibits ferromagnetic alignment.
* Disordered Phase (g > 1): This regime is dominated by the transverse field, and the ground state approaches a product state aligned with the field axis.

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


## Getting Started

### Prerequisites

A Python environment with the following dependencies is required: qiskit, qiskit-ibm-runtime, qiskit-algorithms, numpy, scipy, and matplotlib.

### Running the Diagnostic

To determine the optimal iteration budget for specific hardware or simulator settings, the asymptotic saturation test included in the VQE scripts can be utilized:
1. Set the diagnostic flag to true in the configuration block of the source files.
2. Execute the script to generate an exponential fit of the convergence trajectory to identify the flattening point.

## Key Results Summary

* Optimal Depth: For a system of four qubits, a two-layer hardware-efficient ansatz provides the most effective balance between expressivity and trainability.
* Convergence: Results from the SPSA optimizer reach within approximately 1% error relative to exact diagonalization in ideal simulations when utilizing decaying step sizes.
* Phase Signatures: The framework successfully tracks the collapse of the magnetic structure factor and the peak in entanglement entropy near the critical point.
