import os
import re
import numpy as np
import matplotlib.pyplot as plt

def parse_vector_field(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    match = re.search(r'nonuniform\s+List<vector>\s*\d+\s*\((.*?)\)\s*;', content, re.DOTALL)
    if not match:
        match = re.search(r'\d+\s*\((.*?)\)\s*;', content, re.DOTALL)
    vecs_text = re.findall(r'\((-?\d+\.?\d*e?-?\d*)\s+(-?\d+\.?\d*e?-?\d*)\s+(-?\d+\.?\d*e?-?\d*)\)', match.group(1))
    return np.array([(float(v[0]), float(v[1]), float(v[2])) for v in vecs_text])

def get_latest_time_folder(case_dir):
    time_dirs = []
    for entry in os.listdir(case_dir):
        if os.path.isdir(os.path.join(case_dir, entry)):
            try:
                t = float(entry)
                if t > 0:
                    time_dirs.append((t, entry))
            except ValueError:
                pass
    if not time_dirs:
        return None
    time_dirs.sort()
    return time_dirs[-1][1]

def load_case_fields(case_dir):
    time_folder = get_latest_time_folder(case_dir)
    if not time_folder:
        raise FileNotFoundError(f"No time folder found in {case_dir}")
    
    c_path = os.path.join(case_dir, time_folder, 'C')
    u_path = os.path.join(case_dir, time_folder, 'U')
    
    if not (os.path.exists(c_path) and os.path.exists(u_path)):
        raise FileNotFoundError(f"Missing fields (C or U) in {os.path.join(case_dir, time_folder)}")
        
    C = parse_vector_field(c_path)
    U = parse_vector_field(u_path)
    return C, U

# Case directories
cases = {
    "Grade I (Ks = 0.33 mm)": 'E:/B_ridgi/B_ridgi/Ks_0.33',
    "Grade II (Ks = 0.68 mm)": 'E:/B_ridgi/B_ridgi',
    "Grade III (Ks = 1.90 mm)": 'E:/B_ridgi/B_ridgi/Ks_1.9'
}

fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True, dpi=300)

for idx, (grade, path) in enumerate(cases.items()):
    ax = axes[idx]
    print(f"Loading {grade}...")
    C, U = load_case_fields(path)
    
    x = C[:, 0]
    y = C[:, 1]
    ux = U[:, 0]
    
    # Filter points around the bridge contraction to make the contour focused
    mask = (x >= 0.5) & (x <= 2.5)
    x_f = x[mask]
    y_f = y[mask]
    ux_f = ux[mask]
    
    # Tricontourf to plot the unstructured mesh data
    cntr = ax.tricontourf(x_f, y_f, ux_f, levels=100, cmap='jet')
    
    # Add a grey rectangle representing the solid bridge block
    rect = plt.Rectangle((1.0, 0.075), 0.15, 0.025, facecolor='darkgray', edgecolor='black', hatch='//', zorder=10)
    ax.add_patch(rect)
    
    # Text label for the block
    ax.text(1.075, 0.0825, "Bridge", color='white', weight='bold', fontsize=8, ha='center', va='center', zorder=11)
    
    ax.set_title(grade, fontsize=12, fontweight='bold')
    ax.set_ylabel('y (m)', fontsize=10, fontweight='bold')
    ax.set_xlim(0.5, 2.5)
    ax.set_ylim(0.0, 0.10)
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Add colorbar for each plot
    cbar = fig.colorbar(cntr, ax=ax, orientation='vertical', pad=0.01)
    cbar.set_label('$U_x$ (m/s)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

axes[2].set_xlabel('x (m)', fontsize=10, fontweight='bold')
plt.suptitle('Streamwise Velocity ($U_x$) Contour Comparison', fontsize=16, fontweight='bold', y=0.98)
plt.tight_layout()
output_img = 'E:/B_ridgi/B_ridgi/velocity_contours.png'
plt.savefig(output_img, dpi=300)
print(f"Saved velocity contour comparison to {output_img}")
