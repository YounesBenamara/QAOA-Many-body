# Iteration-Budget Optimization for Stochastic VQE

## 1. Objective

This diagnostic is used to choose a suitable number of iterations for stochastic optimizers such as SPSA in the context of the Variational Quantum Eigensolver (VQE).

Its purpose is to replace an arbitrary iteration count by a data-driven criterion. In practice, it estimates the iteration range beyond which the optimization is no longer improving significantly and becomes dominated by fluctuations. This range is then used to choose both the total iteration budget and the final averaging window.

## 2. Method

Let $E_i(k)$ denote the energy measured at iteration $k$ during the $i$-th optimization run. Because SPSA trajectories are noisy and depend strongly on initialization, a single run is not sufficient to estimate the global convergence trend.

We first average over $M$ independent runs:

$$
\mu(k) = \frac{1}{M} \sum_{i=1}^M E_i(k)
$$

The averaged trajectory is then smoothed with a moving average of width $W$:

$$
\langle E \rangle(k) = \frac{1}{W} \sum_{j=k-W+1}^{k} \mu(j)
$$

To quantify how fast the optimization is still progressing, we define the discrete slope

$$
\Delta(k) = \left| \langle E \rangle(k) - \langle E \rangle(k-1) \right|
$$

In practice, $\Delta(k)$ is often well approximated by the empirical model

$$
\Delta(k) \approx a e^{-bk} + c
$$

where $a$ is the initial transient amplitude, $b$ an effective decay rate, and $c$ the asymptotic floor.

The characteristic iteration $k_{\mathrm{opt}}$ is then defined as the point where the fitted curve reaches a chosen tolerance $\epsilon$:

$$
k_{\mathrm{opt}} = \frac{1}{b} \ln\left(\frac{a}{\epsilon - c}\right)
$$

provided that $c < \epsilon$.

If $c \geq \epsilon$, the threshold is never reached. In that case, we define

$$
k_{\mathrm{opt}} = \frac{1}{b}\ln\left(\frac{a}{c}\right)
$$

which marks the onset of the noise-dominated regime.

## 3. Derived Parameters

Once $k_{\mathrm{opt}}$ is estimated, the production parameters are chosen as follows:

- **Maximum iteration budget**:
  $$
  K_{\max} = \lfloor 1.2 \, k_{\mathrm{opt}} \rfloor
  $$

- **Final averaging window**:
  the Polyak-Ruppert averaging window is chosen within the interval $[k_{\mathrm{opt}}, K_{\max}]$.

- **Effective decay rate $b$**:
  this parameter gives a useful indication of the convergence speed for the chosen optimizer, ansatz, and Hamiltonian.

## 4. Relevance for VQE

This diagnostic is useful in the VQE setting for three main reasons:

1. It avoids choosing the iteration budget arbitrarily.
2. It indicates when further optimization becomes dominated by statistical or hardware noise.
3. It provides a simple way to compare convergence behavior across simulators and hardware backends.

## 5. Interpretation

This method should be viewed as an empirical calibration tool rather than a rigorous convergence test. Its role is to estimate a practical iteration budget for stochastic VQE runs.
