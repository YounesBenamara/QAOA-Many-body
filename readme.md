# Quantum phase transition detection with VQE and QAOA on the transverse-field Ising model

## Overview

This project studies quantum phase transition detection in the one-dimensional **transverse-field Ising model (TFIM)** with **periodic boundary conditions** using three complementary approaches:

1. **Exact diagonalization** as a reference baseline.
2. **VQE (Variational Quantum Eigensolver)** on simulators and IBM quantum hardware.
3. **QAOA (Quantum Approximate Optimization Algorithm)** as a second variational strategy for ground-state approximation and phase-sensitive observables.

The objective is to track how the ground state changes as the transverse field strength is varied, and to determine whether variational quantum algorithms can capture the signatures of the quantum phase transition near the critical point.

---

## Physical model

### Transverse-field Ising Hamiltonian

We consider a ring of $N$ spin-$1/2$ particles with Hamiltonian $H$ modeling the impact of the nearest neighbors interactions and the transverse magnetic field:

$$
H(g) = - \sum_{i=0}^{N-1} Z_i Z_{i+1} - g \sum_{i=0}^{N-1} X_i,
$$

with periodic boundary conditions:

$$
Z_N \equiv Z_0.
$$

Here:

- $Z_i$ and $X_i$ are Pauli operators acting on qubit $i$,
- the term $-Z_i Z_{i+1}$ favors ferromagnetic alignment in the computational basis,
- the transverse-field term $-g X_i$ tends to delocalize the spins and induce quantum fluctuations.

![TFIM](figures/readme_figs/TFIM.jpg)

### Physical regimes

For finite $N$, there is no true singularity, but the model still exhibits a clear crossover which sharpens toward the thermodynamic-limit quantum critical point:

$$
g_c = 1.
$$

The two main regimes are:

- **Ferromagnetic / ordered regime**: $g < 1$
- **Paramagnetic / disordered regime**: $g > 1$

At small $g$, the interaction term dominates and the ground state is close to the ordered state $\ket{\uparrow\uparrow\uparrow...\uparrow}$ or $\ket{\downarrow\downarrow\downarrow...\downarrow}$ with very few spin flips. At large $g$, the transverse field dominates and the ground state approaches an $X$-polarized state $\ket{+++...+}$ where all the spins are aligned with the field.



## Quantities of interest

### 1. Ground-state and first-excited-state energies

The exact solver computes the two lowest eigenvalues:

$$
E_0(g), \qquad E_1(g),
$$

as well as the ground-state eigenvector $|\psi_0(g)\rangle$.

Energy densities are monitored as:

$$
\frac{E_0(g)}{N}, \qquad \frac{E_1(g)}{N}.
$$

The evolution of these quantities with $g$ gives a first signature of the change in physical regime.

### 2. Magnetic structure factor

To quantify spin ordering, the project uses the longitudinal structure factor

$$
S_{zz} = \frac{1}{N^2} \sum_{i,j} \langle Z_i Z_j \rangle.
$$

In the ordered phase, long-range $Z$-correlations are strong, so $S_{zz}$ is large. In the disordered phase, these correlations weaken, so $S_{zz}$ decreases.

### 3. Bipartite von Neumann entropy

For a bipartition of the chain into two subsystems $A$ and $B$, with reduced density matrix

$$
\rho_A = \mathrm{Tr}_B\left(|\psi_0\rangle\langle\psi_0|\right),
$$

we compute the von Neumann entropy

$$
S_{\mathrm{vN}} = -\mathrm{Tr}(\rho_A \log \rho_A).
$$

This quantity is useful because entanglement typically changes significantly around the critical region.

---

## Code architecture

The project is organized around three main layers.

### 1. Exact solver  (`exact_solver.py`)

The exact solver is used to:

- build the TFIM Hamiltonian as a `SparsePauliOp`,
- compute the two lowest eigenvalues using sparse diagonalization,
- obtain the exact ground-state vector,
- derive observables such as structure factor and bipartite entropy.

### 2. VQE

The VQE workflow uses three distinct implementations to progress from conceptual validation to real hardware execution:

- **`vqe_v1.py` (Baseline)**: A streamlined script using Qiskit's native optimizers (COBYLA, SPSA). Designed for clean, predictable baseline comparisons in both ideal and noisy simulated environments.
- **`vqe_v2.py` (Custom Optimizer)**: An advanced implementation featuring a manual SPSA optimization loop. It includes custom decaying learning rates, noise-aware early stopping criteria, fine-tuned parameter initialization, and incremental data saving for deeper algorithmic control.
- **`vqe_v3_hardware.py` (Hardware Execution)**: Production-ready script adapted from `v2`, specifically tailored for deployment on actual IBM Quantum processors (e.g., IBM Marrakech) using Qiskit Runtime Primitives.

The variational circuit is a hardware-efficient ansatz, with a real-valued version based on $R_y$ rotations. This is physically motivated because the TFIM Hamiltonian is real in the computational basis.

### 3. QAOA

The QAOA part of the project is intended to provide a second variational strategy to probe the same Hamiltonian and compare its performance with VQE.

The comparison is relevant because:

- VQE uses a hardware-efficient ansatz with generic variational freedom,
- QAOA uses a structured ansatz directly tied to the Hamiltonian terms.

This may affect convergence, physical interpretability, and robustness near the critical region.

---

## Exact solver: theoretical and implementation details

### Hamiltonian construction

The TFIM Hamiltonian is encoded as a sum of Pauli strings:

```python
sparse_list = []
for i in range(N):
    sparse_list.append(("ZZ", [i, (i + 1) % N], -1.0))
for i in range(N):
    sparse_list.append(("X", [i], -g))
```

This corresponds exactly to

$$
H(g) = - \sum_i Z_i Z_{i+1} - g \sum_i X_i.
$$

### Lowest-energy computation

The project then converts the Hamiltonian to a sparse matrix and computes the two lowest eigenvalues using `eigsh`, which is appropriate for this type of Hermitian sparse problem.

### Structure factor

The implemented structure factor is based on the operator

$$
\sum_{i,j} Z_i Z_j,
$$

normalized by $N^2$.

The code reconstructs this quantity from the exact ground state and uses it to track the decay of ferromagnetic ordering as $g$ increases.

### Bipartite entropy

The project also computes the half-chain von Neumann entropy by tracing out half the qubits. This gives an entanglement-based signature of the phase crossover and is particularly informative near the critical region.

---

## VQE methodology

### Ansatz choice

The VQE workflow uses a hardware-efficient ansatz. In the real-valued setting, the circuit is constructed using $R_y$ layers and entangling CNOTs.

This choice is motivated by two facts:

1. the TFIM Hamiltonian is real,
2. shallower real-valued ansätze are more hardware-friendly than a fully generic rotation structure.

### Why exact diagonalization is necessary

Exact diagonalization provides:

- the true ground-state energy,
- the true first-excited-state energy,
- exact observables,
- a reference to measure the variational error.

Without this baseline, it would be difficult to determine whether a given VQE or QAOA result is physically meaningful.

### SPSA and COBYLA observations so far

The main optimizer-related conclusion reached so far is:

- **COBYLA performs very well in ideal simulation** and reaches near-exact energies.
- An initial SPSA implementation performed poorly, but this was largely due to poor initialization and suboptimal hyperparameters.
- After switching to:
  - a **near-zero initialization** for small $g$ values,
  - **decaying SPSA step sizes**,
  - a **noise-aware plateau criterion**,
  the SPSA results improved substantially.

This is physically consistent with the TFIM at small $g$: the ground state is already close to the ordered $Z$ sector, so initializing near the identity circuit is much more meaningful than starting from a random point in parameter space.

---

## Current results summary

### 1. Exact solver

The exact solver successfully produces:

- $E_0(g)$ and $E_1(g)$,
- energy density plots,
- structure-factor plots,
- entanglement entropy plots.

These results provide the reference signatures of the phase transition and establish the critical region near $g=1$.

### 2. VQE in ideal simulation

The tuned SPSA version in ideal simulation reaches energies extremely close to the exact reference, showing that:

- the ansatz is expressive enough,
- the improved SPSA schedule is now adequate,
- the optimization problem is under control in the noiseless setting.

This validates the VQE workflow itself.

### 3. VQE on fake backend / noisy simulation

On the fake backend, the optimized energy remains above the exact value by a noticeable but much smaller margin than before the SPSA corrections.

This indicates that once the optimizer and initialization are fixed, the remaining discrepancy is mostly due to:

- noise,
- transpilation overhead,
- imperfect hardware-level execution.

This is an important milestone: it means the noisy-backend error can now be interpreted as a realistic degradation rather than as a failure of the variational method itself.

### 4. Real IBM hardware execution

The real-hardware VQE stage is intended to test whether these trends persist on an actual superconducting quantum processor.

Key constraints in this stage are:

- limited runtime budget,
- shot noise,
- gate errors,
- transpilation overhead,
- queue constraints.

The project therefore uses shallow circuits, shot-based estimation, and noise-aware convergence criteria.

---

## Suggested figures to include in this README

### Exact solver section

- **TFIM ring lattice**
- **Ground-state and first-excited-state energy densities vs $g$**
- **Structure factor vs $g$ for several system sizes**
- **Half-chain von Neumann entropy vs $g$ for several system sizes**

### VQE section

- **Ideal VQE convergence curve**
- **Fake-backend VQE convergence curve**
- **Real-hardware VQE convergence curve**
- **Comparison plot: exact vs VQE energy as a function of $g$**
- **Comparison table: COBYLA vs SPSA in ideal simulation**

### QAOA section

- **QAOA convergence curve**
- **QAOA vs VQE comparison for selected values of $g$**
- **Observable comparison near the critical point**

---

## Results section template

## Results

### Exact diagonalization

_Insert here a short discussion of the exact spectra and observables._

**Figure placeholders:**

- `figures/exact_solver_figs/energy_density_N=...png`
- `figures/exact_solver_figs/structure_factor.png`
- `figures/exact_solver_figs/VNeumann_entropy.png`

### VQE — ideal simulation

_Insert here the exact numerical performance of VQE in noiseless mode._

Suggested points to report:

- optimizer used,
- ansatz depth,
- number of parameters,
- final energy,
- relative error with respect to exact diagonalization,
- interpretation of convergence behavior.

### VQE — fake backend / noisy simulation

_Insert here the noisy simulated results._

Suggested discussion points:

- effect of initialization,
- impact of SPSA scheduling,
- remaining gap with the exact energy,
- interpretation in terms of hardware noise and transpilation.

### VQE — real IBM hardware

_Insert here the real QPU results._

Suggested discussion points:

- backend used,
- number of shots,
- transpiled depth and gate counts,
- final energy and uncertainty,
- difference relative to exact and fake-backend results.

### QAOA

_Insert here the QAOA results and their comparison to VQE._

Suggested discussion points:

- QAOA depth $p$,
- final energy,
- optimizer used,
- comparison with exact/VQE,
- whether QAOA captures phase-sensitive observables effectively.

---

## Main conclusions so far

At the current stage, the project supports the following conclusions:

1. The exact TFIM solver provides a reliable baseline for energy and observable benchmarking.
2. The chosen variational ansatz is expressive enough for small-system TFIM ground-state approximation.
3. Optimizer design and initialization matter strongly, especially for SPSA.
4. Once the optimizer is corrected, the remaining gap on noisy backends is largely attributable to hardware-like noise rather than variational failure.
5. This makes the project a suitable framework for comparing VQE and QAOA as tools for quantum phase-transition detection.

---

## Possible next steps

- Extend the VQE study to multiple values of $g$ across the transition region.
- Compare VQE and QAOA not only on energy, but also on observables such as $S_{zz}$ and entanglement proxies.
- Study finite-size effects by repeating the workflow for several values of $N$.
- Improve the real-hardware workflow using backend-aware transpilation and measurement strategies.
- Add a direct section analyzing how well the variational states detect the critical region around $g=1$.

---

## Repository structure

```text
project/
├── src/
│   ├── exact_solver.py
│   ├── vqe_v1.py
│   ├── vqe_v2.py
│   ├── vqe_v3_hardware.py
│   ├── g_comparison.py
│   ├── layers_comparison.py
│   ├── plot.py
│   └── jobs.py
├── figures/
│   ├── exact_solver_figs/
│   ├── vqe_figs/
│   └── qaoa_figs/
└── README.md
```

---

## Dependencies

Typical dependencies used in this project include:

- `qiskit`
- `qiskit-ibm-runtime`
- `numpy`
- `scipy`
- `matplotlib`
- `networkx`
- `qctrlvisualizer`

---

## Notes

This README is structured to serve both as:

- a theoretical introduction to the project,
- and a report scaffold for progressively inserting exact, VQE, and QAOA results.

The next step is to replace the placeholders with the actual generated figures and the numerical results already obtained.

