"""
Simple VQE convergence plot: energy evolution + exact value + std band.
Automatically finds the latest run in runs/vqe_hardware/.
"""

from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt


# Auto-find the latest run directory
RUNS_DIR = Path("runs/vqe_hardware")
run_dir = sorted(RUNS_DIR.iterdir())[-1]
print(f"Plotting: {run_dir}")

arrays = np.load(run_dir / "final_arrays.npz")

cost = arrays["cost_history"]
std = arrays["std_history"]
E0 = float(arrays["exact_energy"][0])
Ef = float(arrays["final_energy"][0])
Ef_se = float(arrays["final_energy_se"][0])

rel_err = abs(Ef - E0) / abs(E0)
x = np.arange(1, len(cost) + 1)

with open(run_dir / "metadata.json", "r") as f:
    md = json.load(f)
phys = md["physics"]
backend = md["backend_name"]

plt.figure(figsize=(9, 5))

plt.errorbar(x, cost, yerr=std, fmt="o-", markersize=4, linewidth=1.5,
             color="#0066CC", ecolor="#0066CC", elinewidth=1, capsize=3,
             label=f"VQE energy (err={rel_err*100:.2f}%)")

plt.axhline(E0, color="#FF3300", linestyle="--", linewidth=1.5,
            label=f"Exact $E_0$ = {E0:.6f}")

plt.xlabel("Iteration")
plt.ylabel(r"Energy $\langle H \rangle$")
plt.title(f"VQE on {backend} — N={phys['N']}, g={phys['g']}, L={phys['n_layers']}")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

out = run_dir / "convergence.png"
plt.savefig(out, dpi=300)
plt.show()
print(f"Saved: {out}")