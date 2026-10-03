import shutil
import subprocess
from pathlib import Path
import csv
import numpy as np


# ============================================================
# SETTINGS
# ============================================================

SIRISMP_EXE = Path("../sirismp")

MESH_FILE = Path("siris-output.off")

# phase matrix INPUT file
PHASE_MATRIX_FILE = "outputS.out"

# ONLY used if EXPERIMENTAL_MFP_ON = 1
MFP_DISTRIBUTION_FILE = "~/dists/cdfconstant04dist.txt"

# fixed wavelength
WAVELEN = 6.283185307179586

# particle sizes to test
MESH_SCALES = [200, 20]

NRAYS = 100000
SEED = 0

DIFFUSE_ON = 1
ALBEDO = 1.0
MEAN_FREE_PATH = 0.946029056271191

EXPERIMENTAL_MFP_ON = 1

BEAM_RADIUS = -1


# ============================================================
# RED LIGHT MATERIALS
# ============================================================

MATERIALS_RED = {
    "Halley-like dust": (1.98, 0.48),
    "Amorphous silicate (forsterite)": (1.677, 0.0044),
    "Water ice": (1.308, 1.43e-8),
    "Carbon dioxide ice": (1.41, 1.05e-6),
    "Iron-rich pyroxene": (1.675, 0.0212),
    "Olivine": (1.743, 0.066),
    "Amorphous carbon": (2.14, 0.805),
    "Pyrrhotite (FeS)": (1.70, 1.86),
    "Iron-rich olivine (fayalite)": (1.85, 0.00077),
    "Cosmic organic refractory": (1.98, 0.2677),
    "Tholin ice": (1.540, 0.0012),
    "Titan tholins": (1.557, 0.0009),
}


# ============================================================
# HELPERS
# ============================================================

def safe_name(name):
    return (
        name.lower()
        .replace(" ", "_")
        .replace("-", "_")
        .replace("(", "")
        .replace(")", "")
    )


def write_input_file(workdir, mesh_scale, n, k):

    input_text = f"""nrays {NRAYS}
max_scattering 200
killswitch_start 70
nbins 180
I_cutoff_limit 0.0000001
seed {SEED}
wavelen {WAVELEN}
mesh_scale {mesh_scale}
mesh {MESH_FILE.name}
force_interaction 1
beam_radius {BEAM_RADIUS}
material1 {n} {k} {DIFFUSE_ON} {ALBEDO} {MEAN_FREE_PATH} {PHASE_MATRIX_FILE} {EXPERIMENTAL_MFP_ON} {MFP_DISTRIBUTION_FILE}
"""

    input_path = workdir / "input.in"
    input_path.write_text(input_text)

    return input_path


def parse_outputS(output_path):

    rows = []

    with open(output_path, "r") as f:
        for line in f:

            line = line.strip()

            if not line:
                continue

            try:
                nums = [float(x) for x in line.split()]

                if np.all(np.isfinite(nums)):
                    rows.append(nums)

            except:
                continue

    data = np.array(rows)

    # column assumptions:
    # x angle S11 S12 ...

    x = data[:, 0]
    s11 = data[:, 1]
    s12 = data[:, 2]

    return x, s11, s12


def find_max_negative_s12_over_s11(output_path):

    x, s11, s12 = parse_outputS(output_path)

    y = -s12 / s11

    valid = np.isfinite(y) & (s11 != 0)

    x = x[valid]
    y = y[valid]

    idx = np.argmax(y)

    return x[idx], y[idx]


# ============================================================
# MAIN
# ============================================================

def main():

    base_dir = Path.cwd()

    runs_dir = base_dir / "batch_runs"
    runs_dir.mkdir(exist_ok=True)

    results = []

    for material_name, (n, k) in MATERIALS_RED.items():

        for mesh_scale in MESH_SCALES:

            run_name = f"{safe_name(material_name)}_{mesh_scale}"

            workdir = runs_dir / run_name
            workdir.mkdir(exist_ok=True)

            shutil.copy(MESH_FILE, workdir / MESH_FILE.name)

            shutil.copy(PHASE_MATRIX_FILE,
                        workdir / PHASE_MATRIX_FILE)

            input_path = write_input_file(
                workdir,
                mesh_scale,
                n,
                k
            )

            print(f"\nRunning {material_name} | mesh_scale={mesh_scale}")

            try:

                subprocess.run(
                    [str(SIRISMP_EXE.resolve()), "input.in"],
                    cwd=workdir,
                    check=True
                )

                output_path = workdir / "outputS.out"

                max_x, max_y = find_max_negative_s12_over_s11(
                    output_path
                )

                print(
                    f"{material_name} | "
                    f"mesh_scale={mesh_scale} | "
                    f"x for max y = {max_x:.6g} | "
                    f"max -S12/S11 = {max_y:.6g}"
                )

                results.append({
                    "material": material_name,
                    "mesh_scale": mesh_scale,
                    "n": n,
                    "k": k,
                    "x_for_max_negative_S12_over_S11": max_x,
                    "max_negative_S12_over_S11": max_y,
                    "run_folder": str(workdir)
                })

            except Exception as e:

                print(
                    f"FAILED: {material_name} | "
                    f"mesh_scale={mesh_scale}"
                )

                print(e)

                results.append({
                    "material": material_name,
                    "mesh_scale": mesh_scale,
                    "n": n,
                    "k": k,
                    "x_for_max_negative_S12_over_S11": None,
                    "max_negative_S12_over_S11": None,
                    "run_folder": str(workdir),
                    "error": str(e)
                })

    if len(results) > 0:

        with open("siris_redlight_results.csv",
                  "w",
                  newline="") as f:

            writer = csv.DictWriter(
                f,
                fieldnames=results[0].keys()
            )

            writer.writeheader()
            writer.writerows(results)

    print("\nDone.")
    print("Saved results to siris_redlight_results.csv")


if __name__ == "__main__":
    main()
