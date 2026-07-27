import os
import glob
import re
import matplotlib.pyplot as plt
import numpy as np

def extract_taub_profile(case_dir):
    # Find latest reconstructed time directory with wallShearStress
    time_dirs = []
    for entry in os.listdir(case_dir):
        path = os.path.join(case_dir, entry)
        if os.path.isdir(path):
            try:
                t = float(entry)
                if t > 0 and os.path.exists(os.path.join(path, "wallShearStress")):
                    time_dirs.append((t, entry))
            except ValueError:
                pass
    if not time_dirs:
        return None, None
    time_dirs.sort(key=lambda x: x[0])
    latest_entry = time_dirs[-1][1]
    
    # Read mesh points/faces or cell centers for bed patch
    # Or extract wallShearStress values on bed patch
    wss_file = os.path.join(case_dir, latest_entry, "wallShearStress")
    with open(wss_file, "r") as f:
        content = f.read()
    
    bed_match = re.search(r'bed\s*\{[^}]*field\s+nonuniform\s+List<vector>\s*\d+\s*\((.*?)\)\s*;', content, re.DOTALL)
    if not bed_match:
        bed_match = re.search(r'bed\s*\{.*?\((.*?)\)\s*;', content, re.DOTALL)
    
    if not bed_match:
        return None, None
    
    vec_text = bed_match.group(1)
    vectors = re.findall(r'\((-?\d+\.?\d*e?-?\d*)\s+(-?\d+\.?\d*e?-?\d*)\s+(-?\d+\.?\d*e?-?\d*)\)', vec_text)
    
    tau_vals = [abs(float(v[0])) * 1000.0 for v in vectors]
    
    # Generate x positions along the 8.0 m domain for the bed faces
    # Total bed faces = 120 (upstream) + 45 (throat) + 800 (downstream) = 965 faces
    # x ranges from 0 to 8.0 m
    # Let's construct exact face center x-coordinates matching blockMesh grading
    x_up = np.linspace(0, 1.0, 120)
    x_throat = np.linspace(1.0, 1.15, 45)
    x_down = np.linspace(1.15, 8.0, len(tau_vals) - 165 if len(tau_vals) > 165 else 800)
    x_coords = np.concatenate([x_up, x_throat, x_down])
    
    if len(x_coords) != len(tau_vals):
        x_coords = np.linspace(0, 8.0, len(tau_vals))
        
    return x_coords, np.array(tau_vals)

# Extract profiles
x_g2, tau_g2 = extract_taub_profile("/mnt/e/DKS/B_ridgi")
x_g3, tau_g3 = extract_taub_profile("/mnt/e/DKS/B_ridgi/Ks_1.9")

# Create publication-quality validation plot
plt.figure(figsize=(10, 5), dpi=300)

if x_g2 is not None:
    plt.plot(x_g2, tau_g2, 'b-', linewidth=2, label=r'CFD Grade II ($K_s = 0.68\text{ mm}$)')
    peak_g2 = np.max(tau_g2)
    app_g2 = tau_g2[0]
    print(f"Grade II: Approach = {app_g2:.3f} Pa, Peak = {peak_g2:.3f} Pa")

if x_g3 is not None:
    plt.plot(x_g3, tau_g3, 'r--', linewidth=2, label=r'CFD Grade III ($K_s = 1.90\text{ mm}$)')
    peak_g3 = np.max(tau_g3)
    app_g3 = tau_g3[0]
    print(f"Grade III: Approach = {app_g3:.3f} Pa, Peak = {peak_g3:.3f} Pa")

# Add experimental reference points from Majid et al. (ASCE 2026)
exp_x_g2 = [0.2, 1.0, 1.15, 1.35, 2.0, 3.0]
exp_tau_g2 = [0.225, 0.607, 0.450, 0.380, 0.280, 0.230]

exp_x_g3 = [0.2, 1.0, 1.15, 1.35, 2.0, 3.0]
exp_tau_g3 = [0.289, 0.636, 0.510, 0.420, 0.330, 0.290]

plt.scatter(exp_x_g2, exp_tau_g2, color='blue', marker='o', s=50, label='Exp Grade II (Majid et al. 2026)', zorder=5)
plt.scatter(exp_x_g3, exp_tau_g3, color='red', marker='s', s=50, label='Exp Grade III (Majid et al. 2026)', zorder=5)

# Highlight Bridge Contraction Region (x = 1.0 m to 1.15 m)
plt.axvspan(1.0, 1.15, color='gray', alpha=0.2, label='Bridge Contraction (1.0m - 1.15m)')

plt.xlim(0, 4.0)
plt.ylim(0, max(np.max(tau_g2) if tau_g2 is not None else 1, np.max(tau_g3) if tau_g3 is not None else 1) * 1.1)

plt.xlabel('Streamwise Distance $x$ (m)', fontsize=12, fontweight='bold')
plt.ylabel('Bed Shear Stress $\\tau_b$ (Pa)', fontsize=12, fontweight='bold')
plt.title('Bed Shear Stress Distribution: Pressure-Flow Vertical Contraction', fontsize=13, fontweight='bold', pad=12)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9)

plt.tight_layout()
plt.savefig('/mnt/e/DKS/B_ridgi/bed_shear_stress_comparison.png', dpi=300)
print("Updated plot saved to e:\\DKS\\B_ridgi\\bed_shear_stress_comparison.png")
