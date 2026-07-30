import os
import re
import json
import matplotlib.pyplot as plt
import numpy as np

# Reference open-channel shear stress values (Pa)
tau_o_vals = {
    "Grade I": 0.137,
    "Grade II": 0.180,
    "Grade III": 0.282
}

def get_bed_taub(case_dir):
    points_file = os.path.join(case_dir, 'constant/polyMesh/points')
    faces_file = os.path.join(case_dir, 'constant/polyMesh/faces')
    boundary_file = os.path.join(case_dir, 'constant/polyMesh/boundary')
    
    if not os.path.exists(points_file):
        return None, None
        
    with open(points_file, 'r') as f:
        pts_text = re.findall(r'\((-?\d+\.?\d*e?-?\d*)\s+(-?\d+\.?\d*e?-?\d*)\s+(-?\d+\.?\d*e?-?\d*)\)', f.read())
    pts = [(float(p[0]), float(p[1]), float(p[2])) for p in pts_text]
    
    with open(boundary_file, 'r') as f:
        bnd_content = f.read()
    bed_bnd = re.search(r'bed\s*\{[^}]*nFaces\s+(\d+)\s*;[^}]*startFace\s+(\d+)\s*;', bnd_content)
    n_faces = int(bed_bnd.group(1))
    start_face = int(bed_bnd.group(2))
    
    with open(faces_file, 'r') as f:
        faces_text = re.findall(r'\d+\((.*?)\)', f.read())
    
    bed_x = []
    for i in range(start_face, start_face + n_faces):
        f_pt_indices = [int(idx) for idx in faces_text[i].split()]
        f_x = sum([pts[idx][0] for idx in f_pt_indices]) / len(f_pt_indices)
        bed_x.append(f_x)
        
    time_dirs = []
    for entry in os.listdir(case_dir):
        path = os.path.join(case_dir, entry)
        if os.path.isdir(path):
            try:
                t = float(entry)
                if t > 0 and os.path.exists(os.path.join(path, 'wallShearStress')):
                    time_dirs.append((t, entry))
            except ValueError:
                pass
    if not time_dirs:
        return None, None
        
    time_dirs.sort(key=lambda item: item[0])
    latest_entry = time_dirs[-1][1]
    
    with open(os.path.join(case_dir, latest_entry, 'wallShearStress'), 'r') as f:
        content = f.read()
    bed_match = re.search(r'bed\s*\{[^}]*field\s+nonuniform\s+List<vector>\s*\d+\s*\((.*?)\)\s*;', content, re.DOTALL)
    if not bed_match:
        bed_match = re.search(r'bed\s*\{.*?\((.*?)\)\s*;', content, re.DOTALL)
    vec_text = bed_match.group(1)
    vectors = re.findall(r'\((-?\d+\.?\d*e?-?\d*)\s+(-?\d+\.?\d*e?-?\d*)\s+(-?\d+\.?\d*e?-?\d*)\)', vec_text)
    tau_vals = [abs(float(v[0])) * 1000.0 for v in vectors]
    
    # Pair and sort by x-coordinate
    pairs = sorted(zip(bed_x, tau_vals), key=lambda item: item[0])
    xs = np.array([p[0] for p in pairs])
    taus = np.array([p[1] for p in pairs])
    return xs, taus

# Paths to data files
cfd_json_path = r"C:\Users\DKS\.gemini\antigravity-ide\brain\d33da93a-2c81-4f0f-b475-1e2aa22fde5d\scratch\cfd_mrssm_data.json"
exp_json_path = r"C:\Users\DKS\.gemini\antigravity-ide\brain\d33da93a-2c81-4f0f-b475-1e2aa22fde5d\scratch\clean_markers.json"

# Load experimental data
try:
    with open(exp_json_path, 'r') as f:
        exp_data = json.load(f)
    exp_x_g1 = [p["x_val"] for p in exp_data["Grade I (Green)"]]
    exp_tau_g1 = [p["y_val"] for p in exp_data["Grade I (Green)"]]
    exp_x_g2 = [p["x_val"] for p in exp_data["Grade II (Blue)"]]
    exp_tau_g2 = [p["y_val"] for p in exp_data["Grade II (Blue)"]]
    exp_x_g3 = [p["x_val"] for p in exp_data["Grade III (Red)"]]
    exp_tau_g3 = [p["y_val"] for p in exp_data["Grade III (Red)"]]
except Exception as e:
    print(f"Failed to load Experimental markers: {e}")
    exp_x_g1, exp_tau_g1 = [], []
    exp_x_g2, exp_tau_g2 = [], []
    exp_x_g3, exp_tau_g3 = [], []

# =========================================================================
# PART 1: Direct Wall Shear Stress Plots
# =========================================================================
print("\n--- Generating Direct Wall Shear Stress Plots ---")
x_wall_g1, tau_wall_g1 = get_bed_taub('E:/B_ridgi/B_ridgi/Ks_0.33')
x_wall_g2, tau_wall_g2 = get_bed_taub('E:/B_ridgi/B_ridgi')
x_wall_g3, tau_wall_g3 = get_bed_taub('E:/B_ridgi/B_ridgi/Ks_1.9')

# Create subplots
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True, dpi=300)
for idx, (x_cfd, tau_cfd, exp_x, exp_tau, color, marker, label, limit_y, ax) in enumerate([
    (x_wall_g1, tau_wall_g1, exp_x_g1, exp_tau_g1, 'green', '^', 'Grade I ($K_s = 0.33\text{ mm}$)', 3.5, axes[0]),
    (x_wall_g2, tau_wall_g2, exp_x_g2, exp_tau_g2, 'blue', 'o', 'Grade II ($K_s = 0.68\text{ mm}$)', 3.5, axes[1]),
    (x_wall_g3, tau_wall_g3, exp_x_g3, exp_tau_g3, 'red', 's', 'Grade III ($K_s = 1.90\text{ mm}$)', 3.5, axes[2])
]):
    if x_cfd is not None:
        x_norm = (x_cfd - 1.0) / 0.1
        tau_norm = tau_cfd / tau_o_vals[list(tau_o_vals.keys())[idx]]
        ax.plot(x_norm, tau_norm, color=color, linestyle='-', linewidth=2, label=f'CFD {label} (Wall)')
        print(f"{list(tau_o_vals.keys())[idx]} (Wall): Peak = {np.max(tau_norm[(x_norm >= 0) & (x_norm <= 5.0)]):.4f} (normalized)")
    if exp_x:
        ax.scatter(exp_x, exp_tau, color=color, marker=marker, s=30, alpha=0.7, label=f'Exp {label}', zorder=5)
    ax.axvspan(0.0, 1.5, color='gray', alpha=0.2, label='Contraction Region' if idx == 0 else "")
    ax.set_xlim(0, 5.0)
    ax.set_ylim(0, 7.5) # Extended limit because wall shear stress goes high
    ax.set_xlabel('Normalized Streamwise Distance $(x - x_0)/H_a$', fontsize=11, fontweight='bold')
    if idx == 0:
        ax.set_ylabel('Normalized Bed Shear Stress $\\tau_b / \\tau_{o}$', fontsize=11, fontweight='bold')
    ax.set_title(label, fontsize=12, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend(loc='upper right', fontsize=8)

plt.suptitle('Normalized Direct Wall Shear Stress Validation (Direct Gradient)', fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()
plt.savefig('E:/B_ridgi/B_ridgi/bed_shear_stress_comparison_subplots_wall.png', dpi=300)

# Generate individual plots
def generate_individual_plot(x_cfd, tau_cfd, exp_x, exp_tau, color, marker, label, title, filename, tau_o, is_wall=True):
    plt.figure(figsize=(8, 5), dpi=300)
    if x_cfd is not None:
        x_norm = (x_cfd - 1.0) / 0.1
        tau_norm = tau_cfd / tau_o
        plt.plot(x_norm, tau_norm, color=color, linestyle='-', linewidth=2, label=f'CFD {label}')
    if exp_x:
        plt.scatter(exp_x, exp_tau, color=color, marker=marker, s=35, alpha=0.7, label=f'Exp {label}', zorder=5)
    plt.axvspan(0.0, 1.5, color='gray', alpha=0.2, label='Contraction Region')
    plt.xlim(0, 5.0)
    plt.ylim(0, 7.5 if is_wall else 3.5)
    plt.xlabel('Normalized Streamwise Distance $(x - x_0)/H_a$', fontsize=12, fontweight='bold')
    plt.ylabel('Normalized Bed Shear Stress $\\tau_b / \\tau_{o}$', fontsize=12, fontweight='bold')
    plt.title(title, fontsize=13, fontweight='bold', pad=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9, fontsize=10)
    plt.tight_layout()
    plt.savefig(f'E:/B_ridgi/B_ridgi/{filename}', dpi=300)
    plt.close()

generate_individual_plot(x_wall_g1, tau_wall_g1, exp_x_g1, exp_tau_g1, 'green', '^', 'Grade I ($K_s = 0.33\text{ mm}$)', 'Normalized Bed Shear Stress: Grade I (Direct Wall)', 'comparison_grade_I_wall.png', tau_o_vals["Grade I"], is_wall=True)
generate_individual_plot(x_wall_g2, tau_wall_g2, exp_x_g2, exp_tau_g2, 'blue', 'o', 'Grade II ($K_s = 0.68\text{ mm}$)', 'Normalized Bed Shear Stress: Grade II (Direct Wall)', 'comparison_grade_II_wall.png', tau_o_vals["Grade II"], is_wall=True)
generate_individual_plot(x_wall_g3, tau_wall_g3, exp_x_g3, exp_tau_g3, 'red', 's', 'Grade III ($K_s = 1.90\text{ mm}$)', 'Normalized Bed Shear Stress: Grade III (Direct Wall)', 'comparison_grade_III_wall.png', tau_o_vals["Grade III"], is_wall=True)


# =========================================================================
# PART 2: MRSSM Evaluated at y = 3.7 mm Plots
# =========================================================================
print("\n--- Generating MRSSM Shear Stress Plots (y = 3.7mm) ---")
try:
    with open(cfd_json_path, 'r') as f:
        cfd_data = json.load(f)
    x_mrs_g1 = np.array(cfd_data["Grade I"]["x"])
    tau_mrs_g1 = np.array(cfd_data["Grade I"]["tau_xy"])
    x_mrs_g2 = np.array(cfd_data["Grade II"]["x"])
    tau_mrs_g2 = np.array(cfd_data["Grade II"]["tau_xy"])
    x_mrs_g3 = np.array(cfd_data["Grade III"]["x"])
    tau_mrs_g3 = np.array(cfd_data["Grade III"]["tau_xy"])
except Exception as e:
    print(f"Failed to load CFD MRSSM data: {e}")
    x_mrs_g1, tau_mrs_g1 = None, None
    x_mrs_g2, tau_mrs_g2 = None, None
    x_mrs_g3, tau_mrs_g3 = None, None

# Create subplots
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True, dpi=300)
for idx, (x_cfd, tau_cfd, exp_x, exp_tau, color, marker, label, limit_y, ax) in enumerate([
    (x_mrs_g1, tau_mrs_g1, exp_x_g1, exp_tau_g1, 'green', '^', 'Grade I ($K_s = 0.33\text{ mm}$)', 3.5, axes[0]),
    (x_mrs_g2, tau_mrs_g2, exp_x_g2, exp_tau_g2, 'blue', 'o', 'Grade II ($K_s = 0.68\text{ mm}$)', 3.5, axes[1]),
    (x_mrs_g3, tau_mrs_g3, exp_x_g3, exp_tau_g3, 'red', 's', 'Grade III ($K_s = 1.90\text{ mm}$)', 3.5, axes[2])
]):
    if x_cfd is not None:
        x_norm = (x_cfd - 1.0) / 0.1
        tau_norm = tau_cfd / tau_o_vals[list(tau_o_vals.keys())[idx]]
        ax.plot(x_norm, tau_norm, color=color, linestyle='-', linewidth=2, label=f'CFD {label} (MRSSM)')
        print(f"{list(tau_o_vals.keys())[idx]} (MRSSM): Peak = {np.max(tau_norm[(x_norm >= 0) & (x_norm <= 5.0)]):.4f} (normalized)")
    if exp_x:
        ax.scatter(exp_x, exp_tau, color=color, marker=marker, s=30, alpha=0.7, label=f'Exp {label}', zorder=5)
    ax.axvspan(0.0, 1.5, color='gray', alpha=0.2, label='Contraction Region' if idx == 0 else "")
    ax.set_xlim(0, 5.0)
    ax.set_ylim(0, 3.5)
    ax.set_xlabel('Normalized Streamwise Distance $(x - x_0)/H_a$', fontsize=11, fontweight='bold')
    if idx == 0:
        ax.set_ylabel('Normalized Bed Shear Stress $\\tau_b / \\tau_{o}$', fontsize=11, fontweight='bold')
    ax.set_title(label, fontsize=12, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend(loc='upper right', fontsize=8)

plt.suptitle('Normalized Bed Shear Stress Validation (MRSSM at y = 3.7 mm)', fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()
plt.savefig('E:/B_ridgi/B_ridgi/bed_shear_stress_comparison_subplots_mrssm.png', dpi=300)

generate_individual_plot(x_mrs_g1, tau_mrs_g1, exp_x_g1, exp_tau_g1, 'green', '^', 'Grade I ($K_s = 0.33\text{ mm}$)', 'Normalized Bed Shear Stress: Grade I (MRSSM at y=3.7mm)', 'comparison_grade_I_mrssm.png', tau_o_vals["Grade I"], is_wall=False)
generate_individual_plot(x_mrs_g2, tau_mrs_g2, exp_x_g2, exp_tau_g2, 'blue', 'o', 'Grade II ($K_s = 0.68\text{ mm}$)', 'Normalized Bed Shear Stress: Grade II (MRSSM at y=3.7mm)', 'comparison_grade_II_mrssm.png', tau_o_vals["Grade II"], is_wall=False)
generate_individual_plot(x_mrs_g3, tau_mrs_g3, exp_x_g3, exp_tau_g3, 'red', 's', 'Grade III ($K_s = 1.90\text{ mm}$)', 'Normalized Bed Shear Stress: Grade III (MRSSM at y=3.7mm)', 'comparison_grade_III_mrssm.png', tau_o_vals["Grade III"], is_wall=False)

print("All plots successfully generated!")
