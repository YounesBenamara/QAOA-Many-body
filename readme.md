# Variational Ising QPT

[span_0](start_span)This project utilizes both classical simulation and real hardware execution to implement a variational quantum-classical framework designed to detect and characterize the quantum phase transition in the one-dimensional transverse-field Ising model[span_0](end_span). [span_1](start_span)By utilizing variational algorithms, the repository tracks the evolution of the ground state across the critical point g = 1, transitioning from a ferromagnetic ordered phase to a paramagnetic disordered phase[span_1](end_span).

## Current Status

This project is currently under development and the following points should be noted:
* [span_2](start_span)The implementation of the Quantum Approximate Optimization Algorithm is currently lacking[span_2](end_span).
* A comprehensive list of formal sources and citations is yet to be incorporated into the documentation.
* The repository layout and file organization are subject to change as the project evolves.
* The figures directory is currently being reorganized and cleaned for better clarity.
* A detailed report documenting the full findings and physical analysis is currently in preparation using LaTeX.

## Project Overview

[span_3](start_span)[span_4](start_span)The repository explores the capabilities of algorithms from the noisy intermediate-scale quantum era to capture many-body physics signatures using three complementary approaches[span_3](end_span)[span_4](end_span):

1. [span_5](start_span)Exact Diagonalization: A reference baseline utilizing SciPy for benchmarking energy and physical observables[span_5](end_span).
2. [span_6](start_span)Variational Quantum Eigensolver: Ground state approximation using a hardware-efficient ansatz with optimized SPSA and COBYLA loops[span_6](end_span).
3. [span_7](start_span)QAOA: A structured variational strategy tied directly to the Hamiltonian terms intended to probe phase-sensitive observables[span_7](end_span).

## Physical Model

[span_8](start_span)The simulation considers a ring of N spins with periodic boundary conditions governed by the following Hamiltonian[span_8](end_span):

$$H(g) = - \sum_{i=0}^{N-1} Z_i Z_{i+1} - g \sum_{i=0}^{N-1} X_i$$

[span_9](start_span)The model exhibits two primary regimes[span_9](end_span):
* [span_10](start_span)Ordered Phase (g < 1): This regime is dominated by nearest-neighbor interactions where the ground state exhibits ferromagnetic alignment[span_10](end_span).
* [span_11](start_span)Disordered Phase (g > 1): This regime is dominated by the transverse field, and the ground state approaches a product state aligned with the field axis[span_11](end_span).

## Repository Structure

The current structure of the project is organized as follows:

* runs/: Contains output data and logs from specific execution environments, such as the fake Marrakesh backend and real IBM hardware.
* src/: Contains the primary source code, including the exact solver, custom optimizers, and various VQE versions.
* src/figures/: Stores generated plots for both exact diagonalization and variational runs.
* [span_12](start_span)src/exact_solver.py: This module builds the Hamiltonian and computes exact energies, magnetic structure factors, and von Neumann entanglement entropy[span_12](end_span).
* [span_13](start_span)src/vqe_v1.py, vqe_v2.py, and vqe_v3_hardware.py: These scripts provide VQE implementations progressing from basic simulators to real quantum hardware deployment[span_13](end_span).
* src/saturation.md: Documentation regarding the asymptotic saturation tests used for optimizer convergence.

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
* [span_14](start_span)Phase Signatures: The framework successfully tracks the collapse of the magnetic structure factor and the peak in entanglement entropy near the critical point[span_14](end_span).
