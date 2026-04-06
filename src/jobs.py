import numpy as np
import json
import datetime as dt
from pathlib import Path
from qiskit_ibm_runtime import QiskitRuntimeService

# Connect to IBM Quantum
service = QiskitRuntimeService(
    channel='ibm_quantum_platform',
    instance='crn:v1:bluemix:public:quantum-computing:us-east:a/a898f4569476426f82148246f0d5a9a5:89ffd5e0-2162-40e9-96c0-5e827879d0ec::'
)

NUM_JOBS_TO_FETCH = 24 
backend_name = "ibm_marrakesh"

print(f"Fetching your last {NUM_JOBS_TO_FETCH} jobs from {backend_name}...")
recent_jobs = service.jobs(limit=NUM_JOBS_TO_FETCH, backend_name=backend_name)
recent_jobs.reverse()

cost_history = []
std_history = []
job_ids = []
theta_history = []

for iteration, job in enumerate(recent_jobs):
    print(f"Processing job {iteration + 1}/{NUM_JOBS_TO_FETCH}: {job.job_id()}")
    if job.status() == "DONE":
        try:
            job_result = job.result()
            pub_center = job_result[2]
            e_center = float(pub_center.data.evs.reshape(-1)[0])
            s_center = float(pub_center.data.stds.reshape(-1)[0])
            
            cost_history.append(e_center)
            std_history.append(s_center)
            job_ids.append(job.job_id())
            
            # Try to get parameters from the inputs if possible
            try:
                pubs = job.inputs.get('pubs', [])
                param_vals = pubs[2][2] # (circuit, observable, parameter_values)
                theta_history.append(np.asarray(param_vals, dtype=float).flatten())
            except Exception:
                pass
                
        except Exception as e:
            print(f"Failed to extract info for {job.job_id()}: {e}")

# Calculate theoretical hyperparameters based on your vqe_v3_hardware.py
ak_history = []
ck_history = []
A_LR = 0.08
A_SHIFT = 10.0
ALPHA = 0.602
C_PERT = 0.10
GAMMA = 0.101

for k in range(len(cost_history)):
    ak = A_LR / ((k + 1 + A_SHIFT) ** ALPHA)
    ck = C_PERT / ((k + 1) ** GAMMA)
    ak_history.append(ak)
    ck_history.append(ck)

# Recreate the run directory structure
PROJECT_ROOT = Path(__file__).resolve().parent.parent
run_name = f"tfim_vqe_hw_N4_g0.2_L1_{backend_name}_RECOVERED"
run_dir = PROJECT_ROOT / "runs" / "vqe_hardware" / run_name
run_dir.mkdir(parents=True, exist_ok=True)

final_energy = cost_history[-1] if cost_history else 0.0
final_energy_se = std_history[-1] if std_history else 0.0
theta_opt = theta_history[-1] if theta_history else np.zeros(8)
exact_E0 = -4.040594

print(f"\nSaving recovered data to: {run_dir}")

# Recreate final_arrays.npz
np.savez(
    run_dir / "final_arrays.npz",
    theta_opt=np.asarray(theta_opt, dtype=float),
    cost_history=np.asarray(cost_history, dtype=float),
    std_history=np.asarray(std_history, dtype=float),
    ak_history=np.asarray(ak_history, dtype=float),
    ck_history=np.asarray(ck_history, dtype=float),
    exact_energy=np.asarray([exact_E0], dtype=float),
    final_energy=np.asarray([final_energy], dtype=float),
    final_energy_se=np.asarray([final_energy_se], dtype=float),
    job_ids=np.asarray(job_ids, dtype=object),
    converged_early=np.asarray([False]),
)

# Recreate summary.json
def to_jsonable(obj):
    if isinstance(obj, np.generic): return obj.item()
    if isinstance(obj, np.ndarray): return obj.tolist()
    return obj
    
summary = {
    "timestamp_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    "run_dir": str(run_dir),
    "backend_name": backend_name,
    "iterations_done": len(cost_history),
    "converged_early": False,
    "exact_energy": exact_E0,
    "vqe_energy": final_energy,
    "vqe_energy_se": final_energy_se,
    "relative_error": abs(final_energy - exact_E0) / abs(exact_E0) if exact_E0 != 0 else 0,
    "job_ids": job_ids,
    "theta_opt": to_jsonable(theta_opt),
    "cost_history": cost_history,
    "std_history": std_history,
    "ak_history": ak_history,
    "ck_history": ck_history,
}

with open(run_dir / "summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2)

print("Recovery finished successfully! You can now use this folder for plotting.")