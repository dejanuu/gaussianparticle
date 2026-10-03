import numpy as np
import subprocess, time, os, shutil, random

# --------------------------
# Write SIRISmp input from template
# --------------------------
def write_input(template_path, output_path, params):
    text = open(template_path).read()
    for key, val in params.items():
        text = text.replace("{" + key + "}", str(val))
    open(output_path, "w").write(text)

# --------------------------
# Run one SIRISmp execution
# --------------------------
def run_sirismp(sirismp_path, input_file, cwd):
    proc = subprocess.Popen(
        [sirismp_path, input_file],
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    out, err = proc.communicate()
    return proc.returncode, out, err

# --------------------------
# Read outputS.out
# --------------------------
def read_outputS(path):
    theta = []
    S11 = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                theta.append(float(parts[0]))
                S11.append(float(parts[1]))
    return np.array(theta), np.array(S11)

# --------------------------
# Loss function comparing target vs model
# --------------------------
def loss_function(out_path, target_theta, target_S11):
    theta, S11 = read_outputS(out_path)

    # interpolate to target angles
    model_interp = np.interp(target_theta, theta, S11)

    return np.mean((model_interp - target_S11) ** 2)

# --------------------------
# Simulated annealing
# --------------------------
def anneal(template_path, sirismp_path, target_theta, target_S11):

    # initial parameter guess
    current = {
        "REF_REAL": 1.5,
        "REF_IMAG": 0.001,
        "DIFFUSE_ON": 1,
        "MFP": 100,
        "ALBEDO": 0.9,
    }

    T = 1.0
    cooling = 0.95
    best = current.copy()
    best_loss = 1e99

    os.makedirs("anneal_runs", exist_ok=True)

    for step in range(50):

        print(f"\n=== Step {step}  T={T:.3f} ===")

        write_input(template_path, "anneal_runs/input.in", current)

        rc, out, err = run_sirismp(sirismp_path, "input.in", "anneal_runs")

        if rc != 0:
            print("SIRISmp ERROR:")
            print(err)
            continue

        E = loss_function("anneal_runs/outputS.out", target_theta, target_S11)
        print("Loss =", E)

        if E < best_loss:
            best_loss = E
            best = current.copy()
            print("*** new best ***")

        # propose new params
        candidate = {
            "REF_REAL": current["REF_REAL"] + np.random.normal(0, 0.05), #1,2
            "REF_IMAG": max(0, current["REF_IMAG"] + np.random.normal(0, 1)),
            #"DIFFUSE_ON": 1,
            "MFP": max(5, current["MFP"] + np.random.normal(0, 5)),
            "ALBEDO": np.clip(current["ALBEDO"] + np.random.normal(0, 0.05), 0, 1),
        }

        write_input(template_path, "anneal_runs/input.in", candidate)
        rc, out, err = run_sirismp(sirismp_path, "input.in", "anneal_runs")

        if rc != 0:
            print("SIRISmp ERROR:")
            print(err)
            continue

        E_new = loss_function("anneal_runs/outputS.out", target_theta, target_S11)

        # Metropolis condition
        if (E_new < E) or (np.random.rand() < np.exp((E - E_new) / T)):
            current = candidate

        T *= cooling

    return best, best_loss


# -----------------------------------------------
# Main execution
# -----------------------------------------------
if __name__ == "__main__":
    # load experimental target
    target = np.load("target.npz")
    target_theta = target["theta"]
    target_S11 = target["S11"]

    best, best_loss = anneal(
        template_path="input_template.in",
        sirismp_path="/home/dnels/master/multi-particle/sirismp",
        target_theta=target_theta,
        target_S11=target_S11,
    )

    print("\n=== BEST PARAMETERS ===")
    print(best)
    print("Best loss =", best_loss)
