# Asymptotic Saturation Test for Stochastic Optimization (VQE)

## 1. Objective of the Method

The Asymptotic Saturation Test is an algorithmic diagnostic procedure designed to rigorously parameterize stochastic optimizers (such as SPSA) within the framework of a Variational Quantum Eigensolver (VQE). 

Its purpose is to replace an empirical approach (arbitrary choice of iteration count) with an analytical stopping criterion. It precisely identifies the moment the algorithm exits its transient descent phase to enter a stationary regime around a local minimum, thereby dictating the optimal iteration budget and the Polyak-Ruppert averaging window.

## 2. Mathematical Formalism

The method relies on analyzing the first derivative of the mean energy trajectory from an ensemble of independent optimizations. 

Let $E_i(k)$ be the energy evaluated during the $i$-th run at iteration $k$. Given the high variance induced by SPSA's stochasticity (and potentially quantum projection noise), analyzing a single trajectory is insufficient.

**Step 1: Ensemble Average**
We calculate the average over a batch of $M$ independent optimizations launched from random initial points to extract the global trend:
$$\mu(k) = \frac{1}{M} \sum_{i=1}^M E_i(k)$$

**Step 2: Temporal Smoothing (Low-pass filter)**
To eliminate the residual high-frequency noise from the ensemble average, we apply a moving average over a window of size $W$:
$$\langle E \rangle(k) = \frac{1}{W} \sum_{j=k-W+1}^k \mu(j)$$

**Step 3: Derivation**
We compute the absolute discrete derivative of this smoothed trajectory:
$$\Delta(k) = \left| \frac{d\langle E \rangle}{dk} \right|$$

**Step 4: Exponential Fit and Analytical Saturation Point**

On an ideal statevector simulator, the derivative $\Delta(k)$ decays exponentially. We fit a parametric model to the derivative curve:
$$\Delta(k) \approx a \cdot e^{-bk} + c$$

where $a$ is the initial descent amplitude, $b$ is the decay rate, and $c$ is the asymptotic noise floor. The fit quality is assessed via the coefficient of determination $R^2$.

The geometric saturation point $k_{\text{flat}}$ is then obtained **analytically** by solving for the iteration where the fitted curve crosses a critical threshold $\epsilon$ (typically $10^{-3}$, corresponding to chemical accuracy):

$$k_{\text{flat}} = \frac{1}{b} \ln\left(\frac{a}{\epsilon - c}\right)$$

This approach is deterministic and immune to the stochastic fluctuations that plague threshold-crossing methods (where a single noisy spike can shift the detection point by tens of iterations).

**Fallback:** If the asymptotic floor $c \geq \epsilon$ (as may occur on a noisy QPU where shot noise dominates), the threshold is never crossed. In this case, we use the point where the exponential component has decayed to the level of the floor itself: $k_{\text{flat}} = \frac{1}{b}\ln(a/c)$, indicating the shot-noise-limited regime.

## 3. Derived Parameters

Once $k_{\text{flat}}$ is identified via this diagnostic, the production parameters are set as follows:

* **Maximum Iteration Budget ($K_{\max}$)**: Defined with a safety margin beyond the saturation point, typically $K_{\max} = \lfloor 1.2 \times k_{\text{flat}} \rfloor$. This guarantees that all runs reach the stationary regime without wasting computational resources.
* **Extraction Window (Polyak-Ruppert)**: The size of the final averaging window, used to evaluate the true minimum of a run, is set to strictly cover the range $[k_{\text{flat}}, K_{\max}]$. Taking a wider window would include the descent phase and introduce a positive bias.
* **Decay Rate $b$**: This parameter itself is physically meaningful — it characterizes how quickly Qiskit's SPSA auto-calibration converges on the given ansatz/Hamiltonian combination. A larger $b$ indicates a well-conditioned energy landscape.

## 4. Algorithmic and Physical Relevance

This method is particularly relevant in the quantum computing context for several reasons:

1.  **Prevention of Arbitrary Truncation:** The energy landscapes of hardware-efficient ansätze (HEA) are highly non-convex. The "burn-in" time varies drastically depending on the initial position. This test ensures the algorithm is not artificially interrupted mid-descent, which would heavily bias the evaluation of the minima.
2.  **Identification of the Shot-Noise Limit:** On a QPU or a noisy simulator, the test physically identifies the precision limit. The derivative $\Delta(k)$ will not tend to zero but will hit a floor dominated by the standard error of the measurement (proportional to $1/\sqrt{N_{\text{shots}}}$). Continuing optimization beyond this point is analytically useless, as SPSA gradient fluctuations become indistinguishable from hardware noise.
3.  **Simulator Validation:** Analyzing the tail of $\Delta(k)$ immediately discriminates the type of environment. A continuous exponential decay characterizes a statevector (ideal) simulation, whereas horizontal flattening characterizes an environment subjected to projective or hardware noise.