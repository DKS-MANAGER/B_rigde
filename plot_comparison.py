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

# ============================================================
# Extract spatial profiles for all 3 cases
# ============================================================
x_g1, tau_g1 = get_bed_taub('/mnt/e/DKS/B_ridgi/Ks_0.33')
x_g2, tau_g2 = get_bed_taub('/mnt/e/DKS/B_ridgi')
x_g3, tau_g3 = get_bed_taub('/mnt/e/DKS/B_ridgi/Ks_1.9')

# ============================================================
# Geometry reference:
#   Bridge starts at x_mesh = 1.0 m, ends at x_mesh = 1.15 m
#   Ha = 0.10 m (approach flow depth)
#   H  = 0.025 m (contraction height, H/Ha = 0.25)
#   L  = 0.15 m  (contraction length, L/Ha = 1.5)
#
# Paper coordinate: x/Ha measured from START of contraction
#   x/Ha = (x_mesh - 1.0) / 0.10
#   x_mesh = x/Ha * 0.10 + 1.0
# ============================================================
Ha = 0.10  # m

# ============================================================
# Digitized experimental data from Majid et al. (2026)
# Fig. 9(d): L/Ha = 1.5, H/Ha = 0.25
# x-axis: x/Ha (0 to 5), y-axis: tau_uvm / (tau_uvm)_o
#
# The paper presents NORMALIZED shear stress ratios.
# Our CFD outputs ABSOLUTE shear stress in Pa.
# To compare, we normalize the CFD data by the upstream 
# approach value (mean of tau for x < 0.8 m).
# ============================================================

# Digitized from Fig. 9(d) of Majid et al. (2026)
# L1H2V2D2: V=0.26 m/s, D=0.33 mm (Grade I, ks+ ~ 3.56, hydraulically smooth)
exp_xHa_g1 = np.array([-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
exp_tau_norm_g1 = np.array([1.00, 1.00, 1.00, 1.02, 1.10, 1.50, 1.65, 1.45, 1.30, 1.60, 1.85, 1.95, 1.80, 1.55, 1.35])

# L1H2V2D3: V=0.26 m/s, D=0.68 mm (Grade II, ks+ ~ 2.43, hydraulically smooth)
exp_xHa_g2 = np.array([-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
exp_tau_norm_g2 = np.array([1.00, 1.00, 1.00, 1.03, 1.15, 1.60, 1.80, 1.55, 1.35, 1.70, 2.00, 2.10, 1.95, 1.60, 1.40])

# L1H2V2D4: V=0.26 m/s, D=1.90 mm (Grade III, ks+ ~ 1.48, hydraulically transitional)
exp_xHa_g3 = np.array([-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
exp_tau_norm_g3 = np.array([1.00, 1.00, 1.00, 1.05, 1.25, 1.85, 2.15, 1.80, 1.50, 1.80, 2.15, 2.25, 2.10, 1.70, 1.45])

# Convert experimental x/Ha to mesh x (m)
exp_x_g1 = exp_xHa_g1 * Ha + 1.0
exp_x_g2 = exp_xHa_g2 * Ha + 1.0
exp_x_g3 = exp_xHa_g3 * Ha + 1.0

# ============================================================
# PLOT 1: Absolute Bed Shear Stress (Pa) vs x (m)
# ============================================================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), dpi=300)

# --- Top panel: Absolute shear stress ---
if x_g1 is not None:
    ax1.plot(x_g1, tau_g1, 'g:', linewidth=2, label=r'CFD Grade I ($K_s = 0.33$ mm)')
    app_g1 = np.mean(tau_g1[x_g1 < 0.8])
    peak_g1 = np.max(tau_g1[(x_g1 >= 0.95) & (x_g1 <= 1.20)])
    print(f"Grade I  (0.33mm): Approach = {app_g1:.4f} Pa | Contraction Peak = {peak_g1:.4f} Pa | Ratio = {peak_g1/app_g1:.2f}")

if x_g2 is not None:
    ax1.plot(x_g2, tau_g2, 'b-', linewidth=2, label=r'CFD Grade II ($K_s = 0.68$ mm)')
    app_g2 = np.mean(tau_g2[x_g2 < 0.8])
    peak_g2 = np.max(tau_g2[(x_g2 >= 0.95) & (x_g2 <= 1.20)])
    print(f"Grade II (0.68mm): Approach = {app_g2:.4f} Pa | Contraction Peak = {peak_g2:.4f} Pa | Ratio = {peak_g2/app_g2:.2f}")

if x_g3 is not None:
    ax1.plot(x_g3, tau_g3, 'r--', linewidth=2, label=r'CFD Grade III ($K_s = 1.90$ mm)')
    app_g3 = np.mean(tau_g3[x_g3 < 0.8])
    peak_g3 = np.max(tau_g3[(x_g3 >= 0.95) & (x_g3 <= 1.20)])
    print(f"Grade III(1.90mm): Approach = {app_g3:.4f} Pa | Contraction Peak = {peak_g3:.4f} Pa | Ratio = {peak_g3/app_g3:.2f}")

ax1.axvspan(1.0, 1.15, color='gray', alpha=0.2, label='Bridge Contraction')
ax1.set_xlim(0, 3.5)
ax1.set_ylim(0, None)
ax1.set_xlabel('Streamwise Distance $x$ (m)', fontsize=12, fontweight='bold')
ax1.set_ylabel(r'Bed Shear Stress $\tau_b$ (Pa)', fontsize=12, fontweight='bold')
ax1.set_title('(a) Absolute Bed Shear Stress Distribution', fontsize=13, fontweight='bold', pad=10)
ax1.grid(True, linestyle='--', alpha=0.5)
ax1.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9, fontsize=9)

# --- Bottom panel: Normalized shear stress (CFD + Exp) ---
if x_g1 is not None and app_g1 > 0:
    ax2.plot(x_g1, tau_g1 / app_g1, 'g:', linewidth=2, label=r'CFD Grade I ($K_s = 0.33$ mm)')
if x_g2 is not None and app_g2 > 0:
    ax2.plot(x_g2, tau_g2 / app_g2, 'b-', linewidth=2, label=r'CFD Grade II ($K_s = 0.68$ mm)')
if x_g3 is not None and app_g3 > 0:
    ax2.plot(x_g3, tau_g3 / app_g3, 'r--', linewidth=2, label=r'CFD Grade III ($K_s = 1.90$ mm)')

# Plot experimental data (normalized)
ax2.scatter(exp_x_g1, exp_tau_norm_g1, color='green', marker='^', s=50, 
            label='Exp Grade I (Majid et al. 2026)', zorder=5, edgecolors='darkgreen', linewidths=0.5)
ax2.scatter(exp_x_g2, exp_tau_norm_g2, color='blue', marker='o', s=50, 
            label='Exp Grade II (Majid et al. 2026)', zorder=5, edgecolors='darkblue', linewidths=0.5)
ax2.scatter(exp_x_g3, exp_tau_norm_g3, color='red', marker='s', s=50, 
            label='Exp Grade III (Majid et al. 2026)', zorder=5, edgecolors='darkred', linewidths=0.5)

ax2.axvspan(1.0, 1.15, color='gray', alpha=0.2, label='Bridge Contraction')
ax2.axhline(y=1.0, color='black', linestyle='-.', linewidth=0.8, alpha=0.5)
ax2.set_xlim(0, 3.5)
ax2.set_ylim(0, 4.0)
ax2.set_xlabel('Streamwise Distance $x$ (m)', fontsize=12, fontweight='bold')
ax2.set_ylabel(r'$\tau_b / \tau_{b,o}$ (Normalized)', fontsize=12, fontweight='bold')
ax2.set_title(r'(b) Normalized Bed Shear Stress — CFD vs Experiment [Fig. 9(d), Majid et al. 2026]', fontsize=13, fontweight='bold', pad=10)
ax2.grid(True, linestyle='--', alpha=0.5)
ax2.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9, fontsize=8, ncol=2)

plt.tight_layout()
plt.savefig('/mnt/e/DKS/B_ridgi/bed_shear_stress_comparison.png', dpi=300)
print("\nUpdated dual-panel plot saved to: bed_shear_stress_comparison.png")
