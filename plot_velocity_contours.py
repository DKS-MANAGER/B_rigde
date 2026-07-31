import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

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

# Hydraulic constants
Ha = 0.10       # Approach depth (m)
H = 0.025       # Contraction block height (m)
Va = 0.26       # Approach velocity (m/s)

# Case directories and outputs
cases = {
    "Grade I (Ks = 0.33 mm)": ('E:/B_ridgi/B_ridgi/Ks_0.33', 'velocity_contour_grade_I.png'),
    "Grade II (Ks = 0.68 mm)": ('E:/B_ridgi/B_ridgi', 'velocity_contour_grade_II.png'),
    "Grade III (Ks = 1.90 mm)": ('E:/B_ridgi/B_ridgi/Ks_1.9', 'velocity_contour_grade_III.png')
}

for grade, (path, filename) in cases.items():
    print(f"Plotting styled contour for {grade}...")
    C, U = load_case_fields(path)
    
    # Extract x and y coordinates, and Ux velocity
    x = C[:, 0]
    y = C[:, 1]
    ux = U[:, 0]
    
    # Normalize coordinates
    x_norm = (x - 1.0) / Ha
    y_norm = y / H
    ux_norm = ux / Va
    
    # Filter data to the plot region x_norm in [-0.5, 5.0]
    mask = (x_norm >= -0.5) & (x_norm <= 5.0)
    x_f = x_norm[mask]
    y_f = y_norm[mask]
    ux_f = ux_norm[mask]
    
    # Create triangulation for plotting
    triang = tri.Triangulation(x_f, y_f)
    
    # Set up the figure
    fig, ax = plt.subplots(figsize=(10, 4), dpi=300)
    
    # Plot filled contours (velocity field)
    # Use levels from -0.3 to 1.7 to match the reference colorbar
    levels = np.linspace(-0.3, 1.7, 101)
    cntr_f = ax.tricontourf(triang, ux_f, levels=levels, cmap='jet', extend='both')
    
    # Plot line contours (isolines) with custom steps
    line_levels = [-0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7]
    cntr_l = ax.tricontour(triang, ux_f, levels=line_levels, colors='gray', linewidths=0.5)
    
    # Add inline labels on the contour lines
    ax.clabel(cntr_l, inline=True, fontsize=8, fmt='%.1f', colors='black')
    
    # Draw top boundary/ceiling as grey rectangles (where there is no bridge block)
    # The bridge block itself spans x_norm in [0, 1.5] and y_norm in [3, 4]
    # Ceiling upstream: x_norm in [-0.5, 0], y_norm in [3, 4]
    # Ceiling downstream: x_norm in [1.5, 5.0], y_norm in [3, 4]
    ceil_upstream = plt.Rectangle((-0.5, 3.0), 0.5, 1.0, facecolor='grey', edgecolor='none', zorder=8)
    ceil_downstream = plt.Rectangle((1.5, 3.0), 3.5, 1.0, facecolor='grey', edgecolor='none', zorder=8)
    ax.add_patch(ceil_upstream)
    ax.add_patch(ceil_downstream)
    
    # Draw the solid bridge block as a black rectangle
    bridge_block = plt.Rectangle((0.0, 3.0), 1.5, 1.0, facecolor='black', edgecolor='black', zorder=10)
    ax.add_patch(bridge_block)
    
    # Add labels and formatting
    ax.set_xlabel('$x / H_a$', fontsize=11, fontweight='bold')
    ax.set_ylabel('$y / H$', fontsize=11, fontweight='bold')
    ax.set_xlim(-0.5, 5.0)
    ax.set_ylim(0.0, 4.0)
    
    # Add minor/major ticks to match style
    ax.set_xticks(np.arange(-0.5, 5.1, 0.5))
    ax.set_yticks(np.arange(0, 5, 1))
    ax.tick_params(direction='in', top=True, right=True)
    
    # Add colorbar
    cbar = fig.colorbar(cntr_f, ax=ax, ticks=[-0.3, 0.1, 0.5, 0.9, 1.3, 1.7], pad=0.015, aspect=15)
    cbar.set_label('$\overline{u}/V$', fontsize=11, fontweight='bold')
    cbar.ax.tick_params(labelsize=9)
    
    plt.title(f"Streamwise Velocity Contour - {grade}", fontsize=12, fontweight='bold', pad=10)
    plt.tight_layout()
    
    output_path = f'E:/B_ridgi/B_ridgi/{filename}'
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Successfully saved styled contour to {output_path}")
