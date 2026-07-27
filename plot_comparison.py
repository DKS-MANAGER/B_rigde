import os
import re
import matplotlib.pyplot as plt
import numpy as np

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

# Extract spatial profiles for all 3 cases
x_g1, tau_g1 = get_bed_taub('/mnt/e/DKS/B_ridgi/Ks_0.33')
x_g2, tau_g2 = get_bed_taub('/mnt/e/DKS/B_ridgi')
x_g3, tau_g3 = get_bed_taub('/mnt/e/DKS/B_ridgi/Ks_1.9')

# Create publication-quality validation plot
plt.figure(figsize=(10, 5), dpi=300)

if x_g1 is not None:
    plt.plot(x_g1, tau_g1, 'g:', linewidth=2, label=r'CFD Grade I ($K_s = 0.33\text{ mm}$)')
    app_g1 = np.mean(tau_g1[x_g1 < 0.8])
    peak_g1 = np.max(tau_g1[(x_g1 >= 0.95) & (x_g1 <= 1.20)])
    print(f"Grade I (0.33mm): Approach = {app_g1:.4f} Pa | Contraction Peak = {peak_g1:.4f} Pa")

if x_g2 is not None:
    plt.plot(x_g2, tau_g2, 'b-', linewidth=2, label=r'CFD Grade II ($K_s = 0.68\text{ mm}$)')
    app_g2 = np.mean(tau_g2[x_g2 < 0.8])
    peak_g2 = np.max(tau_g2[(x_g2 >= 0.95) & (x_g2 <= 1.20)])
    print(f"Grade II (0.68mm): Approach = {app_g2:.4f} Pa | Contraction Peak = {peak_g2:.4f} Pa")

if x_g3 is not None:
    plt.plot(x_g3, tau_g3, 'r--', linewidth=2, label=r'CFD Grade III ($K_s = 1.90\text{ mm}$)')
    app_g3 = np.mean(tau_g3[x_g3 < 0.8])
    peak_g3 = np.max(tau_g3[(x_g3 >= 0.95) & (x_g3 <= 1.20)])
    print(f"Grade III (1.90mm): Approach = {app_g3:.4f} Pa | Contraction Peak = {peak_g3:.4f} Pa")

# Add experimental reference points from Majid et al. (ASCE 2026)
exp_x_g1 = [0.2, 1.0, 1.15, 1.35, 2.0, 3.0]
exp_tau_g1 = [0.144, 0.461, 0.350, 0.280, 0.210, 0.170]

exp_x_g2 = [0.2, 1.0, 1.15, 1.35, 2.0, 3.0]
exp_tau_g2 = [0.225, 0.607, 0.450, 0.380, 0.280, 0.230]

exp_x_g3 = [0.2, 1.0, 1.15, 1.35, 2.0, 3.0]
exp_tau_g3 = [0.289, 0.636, 0.510, 0.420, 0.330, 0.290]

plt.scatter(exp_x_g1, exp_tau_g1, color='green', marker='^', s=45, label='Exp Grade I (Majid et al. 2026)', zorder=5)
plt.scatter(exp_x_g2, exp_tau_g2, color='blue', marker='o', s=45, label='Exp Grade II (Majid et al. 2026)', zorder=5)
plt.scatter(exp_x_g3, exp_tau_g3, color='red', marker='s', s=45, label='Exp Grade III (Majid et al. 2026)', zorder=5)

# Highlight Bridge Contraction Region (x = 1.0 m to 1.15 m)
plt.axvspan(1.0, 1.15, color='gray', alpha=0.2, label='Bridge Contraction (1.0m - 1.15m)')

plt.xlim(0, 4.0)
plt.ylim(0, 1.0)

plt.xlabel('Streamwise Distance $x$ (m)', fontsize=12, fontweight='bold')
plt.ylabel('Bed Shear Stress $\\tau_b$ (Pa)', fontsize=12, fontweight='bold')
plt.title('Bed Shear Stress Distribution: Pressure-Flow Vertical Contraction', fontsize=13, fontweight='bold', pad=12)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9, fontsize=9)

plt.tight_layout()
plt.savefig('/mnt/e/DKS/B_ridgi/bed_shear_stress_comparison.png', dpi=300)
print("Updated plot saved to e:\\DKS\\B_ridgi\\bed_shear_stress_comparison.png")
