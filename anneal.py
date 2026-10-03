import numpy as np, subprocess, os, shutil

def run_siris(params):
    """Write input file, run siris1p, parse output, return chi-squared."""
    write_input(params)  # update single-particle-input.in
    subprocess.run(["./siris1p", "single-particle-input.in"], ...)
    return parse_output_chi2("output.out", observed_data)

T = T_start
params = params_init
E = run_siris(params)

for step in range(n_steps):
    # Propose a random perturbation
    params_new = params + np.random.randn(n_params) * step_size
    
    E_new = run_siris(params_new)
    dE = E_new - E
    
    # Metropolis acceptance
    if dE < 0 or np.random.rand() < np.exp(-dE / T):
        params = params_new
        E = E_new
    
    T *= cooling_rate  # e.g. 0.995
