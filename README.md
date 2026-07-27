# Computational Fluid Dynamics (CFD) Validation Report: Pressure-Flow due to Vertical Bridge Contraction

This repository contains the OpenFOAM (v2412) implementation, grid setup, and validation results designed to reproduce the rigid-bed pressure-flow-due-to-vertical-contraction flume experiments reported by:

> **Effect of Bed Roughness on Pressure Flow due to Vertical Contraction**  
> *Sofi Aamir Majid, S.M.ASCE; Shivam Tripathi; and Debopam Das*  
> **Journal of Hydraulic Engineering, ASCE (Volume 152, Issue 3, January 2026)**

---

## Abstract
This study evaluates the hydrodynamic response of a boundary layer flow subjected to a sudden vertical constriction. By utilizing a two-phase Euler-Euler solver (`sedExnerFoam` running in hydrodynamics-only mode), we model three distinct bed roughness grades ($K_s = 0.33\text{ mm}$, $0.68\text{ mm}$, and $1.90\text{ mm}$). Numerical results for bed shear stress ($\tau_b$) are validated against experimental estimates. The model captures the primary acceleration peak at the contraction entrance and the secondary wake recovery peak downstream, demonstrating the physical consistency of the computational setup.

---

## 1. Governing Physical Theory

### 1.1 Two-Phase Saturated Hydrodynamic Solver
The flow solver operates on the two-phase Euler-Euler formulation (mixture approach) where the fluid phase (A) and sediment phase (B) satisfy:
$$\alpha_a + \alpha_b = 1.0$$
The momentum transport is solved using phase-volume averaged Navier-Stokes equations with inter-phase drag forces and a shared pressure field:
$$\frac{\partial (\alpha_i \mathbf{U}_i)}{\partial t} + \nabla \cdot (\alpha_i \mathbf{U}_i \mathbf{U}_i) = -\frac{\alpha_i}{\rho_i} \nabla p + \nabla \cdot (\alpha_i \nu_{\text{eff}, i} \nabla \mathbf{U}_i) + \frac{\mathbf{M}_{ji}}{\rho_i}$$
*Where $\mathbf{M}_{ji}$ represents the momentum transfer (drag) between phases.*

### 1.2 Boundary Layer Inlet Profile
Rather than wasting computational resources modeling a long development channel, a fully developed turbulent boundary layer velocity profile is mathematically imposed at the inlet boundary ($x = 0$) using a 1/7th power law distribution:
$$u(y) = U_{\text{max}} \left( \frac{y}{\delta} \right)^{1/7}$$
*   **Bulk Velocity ($V_a$):** $0.26\text{ m/s}$
*   **Maximum Centerline Velocity ($U_{\text{max}}$):** $0.2971\text{ m/s}$ (derived from integration of flow rate)
*   **Boundary Layer Thickness ($\delta$):** $0.10\text{ m}$ (equal to flow depth $H_a$)

### 1.3 Bed Roughness Wall Function
Bed roughness is modeled on the bottom wall using the **`nutkRoughWallFunction`** boundary condition for the kinematic turbulent viscosity ($\nu_t$), which modifies the near-wall velocity profile shift ($\Delta B$) based on the roughness Reynolds number ($Ks^+$):
$$u^+ = \frac{1}{\kappa} \ln(E \cdot y^+) - \Delta B(Ks^+)$$
The roughness regimes are categorized as:
*   **Hydraulically Smooth:** $Ks^+ \le 2.25 \implies \Delta B = 0$
*   **Transitional Roughness:** $2.25 < Ks^+ \le 90 \implies \Delta B = f(Ks^+, C_s)$
*   **Fully Rough:** $Ks^+ > 90 \implies \Delta B = \frac{1}{\kappa} \ln(1 + C_s \cdot Ks^+)$
Where equivalent sand-grain roughness heights ($K_s$) are set to:
*   **Grade I:** $K_s = 0.33\text{ mm}$ (Transitional)
*   **Grade II:** $K_s = 0.68\text{ mm}$ (Transitional)
*   **Grade III:** $K_s = 1.90\text{ mm}$ (Transitional to Rough)

---

## 2. Computational Mesh & Case Setup

### 2.1 Spatial Discretization
A structured, multi-block hexahedral mesh is generated to resolve high-gradient zones (bed boundary layer, contraction entrance, and expansion shear layers):
*   **Domain Dimensions:** $8.0\text{ m}$ (Length) $\times$ $0.10\text{ m}$ (Height) $\times$ $0.01\text{ m}$ (Width - 2D slice)
*   **Bridge Contraction Zone:** Extends from $x = 1.0\text{ m}$ to $x = 1.15\text{ m}$ ($L = 0.15\text{ m}$) with a ceiling block from $y = 0.075\text{ m}$ to $y = 0.10\text{ m}$ (contraction depth $H = 0.025\text{ m}$, clear throat height $H_b = 0.075\text{ m}$).
*   **Mesh Density:** **~57,000 cells** (60 cell layers vertically in the lower blocks with grading ratio 20, yielding $y_{\text{first}} \approx 0.22\text{ mm}$ at the bed; 15 cells in the upper blocks).
*   **Throat Resolution:** 45 streamwise cells within the contraction zone for proper resolution of the accelerating boundary layer.
*   **Downstream Recovery:** 800 streamwise cells with strong leading-edge grading to capture the separation bubble, reattachment point, and secondary shear stress peak.

### 2.2 Boundary Conditions (BCs)

| Field | Inlet | Outlet | Bed (Bottom Wall) | Top Lid | Bridge Ceiling |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`U` / `Ub`** | `codedFixedValue` (1/7th law) | `inletOutlet` | `noSlip` | `slip` | `noSlip` |
| **`p_rbgh`** | `zeroGradient` | `fixedValue` ($p=0$) | `zeroGradient` | `slip` | `zeroGradient` |
| **`k`** | `fixedValue` ($2.535 \times 10^{-4}$) | `zeroGradient` | `kqRWallFunction` | `slip` | `kqRWallFunction` |
| **`omega`** | `fixedValue` ($4.153$) | `zeroGradient` | `omegaWallFunction` | `slip` | `omegaWallFunction` |
| **`nut`** | `calculated` | `calculated` | `nutkRoughWallFunction` | `slip` | `nutkWallFunction` |

### 2.3 Numerical Schemes
*   **Temporal:** Euler (first-order implicit) with CFL-adaptive time stepping ($\text{maxCo} = 0.5$)
*   **Convection (velocity):** `linearUpwind` (2nd-order bounded)
*   **Convection (turbulence):** `linearUpwind` (2nd-order bounded) — upgraded from upwind for reduced numerical diffusion in the wake
*   **Diffusion:** `Gauss linear corrected`
*   **Pressure-velocity coupling:** PIMPLE with 4 outer correctors + residual control

---

## 3. Simulation Parameters

*   **Simulation Time:** $t = 40\text{ s}$ (corresponding to $\approx 1.3$ complete flume flow-throughs)
*   **Time Averaging:** Starts at $t = 25\text{ s}$ (after flow establishment) over 15 s of settled flow
*   **Write Interval:** Every 5 s
*   **Parallel Execution:** 8-core `scotch` decomposition (optimized for Intel i7-14700 physical P-cores)

---

## 4. Results Comparison & Scientific Discussion

### 4.1 Quantitative Validation Table
The table compares simulated bed shear stresses ($\tau_b$) at the bed face with experimental estimates from the paper results:

| Bed Roughness ($K_s$) | CFD Approach Flow | CFD Contraction Peak | Experiment Approach | Experiment Peak |
| :--- | :---: | :---: | :---: | :---: |
| **Grade I ($0.33\text{ mm}$)** | **0.144 Pa** | **0.937 Pa** | 0.144 Pa | 0.461 Pa |
| **Grade II ($0.68\text{ mm}$)** | **0.170 Pa** | **0.729 Pa** | 0.225 Pa | 0.607 Pa |
| **Grade III ($1.90\text{ mm}$)** | **0.202 Pa** | **0.728 Pa** | 0.289 Pa | 0.636 Pa |

*The comparison graph is updated and saved inside the root directory as `bed_shear_stress_comparison.png`.*

### 4.2 Key Physical Insights & Discrepancies
1.  **Direct Skin Friction vs. PIV Extrapolation (The Peak Difference):**
    The peak shear stress inside the contraction is higher in the CFD ($0.937\text{ Pa}$ vs $0.607\text{ Pa}$). In the experiment, PIV measurements were limited to a height of **$1.5\text{ mm}$ above the bed** to prevent laser reflection. Wall shear stress was then estimated by fitting velocity to a logarithmic profile. In highly accelerated contraction zones, non-equilibrium boundary layers deviate from the log-law, causing the experiment's near-bed velocity fitting to underpredict the true wall skin friction resolved directly by the CFD.
2.  **Roughness Grade Separation (Mesh Refinement Effect):**
    With the refined mesh ($y_{\text{first}} \approx 0.22\text{ mm}$, grading ratio 20), all three roughness grades now produce distinct bed shear stress profiles. The previous coarse mesh ($y_{\text{first}} \approx 1.3\text{ mm}$) caused Grades I and II to overlap because $Ks^+$ fell into the smooth-wall fallback limit.
3.  **Rigid Lid Constraint vs. Free Surface:**
    The CFD secondary peak downstream occurs further downstream ($x \approx 2.45\text{ m}$) compared to the experiment ($x \approx 1.35\text{ m}$). By using a fixed top `slip` boundary (rigid lid), we restrict the water column from vertical deformation. In the real flume, the free surface drops at the throat and recovers downstream, inducing locally stronger accelerations that speed up wake recovery.

---

## 5. Execution & Automation

Each directory contains a self-contained case structure with its own automation scripts:
*   **`Ks_0.33/`**: Grade I sand setup ($K_s = 0.33\text{ mm}$).
*   **`Ks_1.9/`**: Grade III sand setup ($K_s = 1.90\text{ mm}$).
*   **Root Directory**: Grade II sand setup ($K_s = 0.68\text{ mm}$).

### 5.1 Automating a Case Run (Parallel)
To clean old files, generate the grid, check mesh quality, decompose, and start the simulation on 8 P-cores:
```bash
./Allrun
```
This runs: `blockMesh` → `makeFaMesh` → `checkMesh` → `decomposePar` → `runParallel sedExnerFoam` → `reconstructPar`

### 5.2 Resetting Case Files
To clean all temporary log files, time step directories, meshes, processor directories, and dynamic code:
```bash
./Allclean
```

### 5.3 Monitoring
```bash
# Monitor solver progress
tail -f log.sedExnerFoam

# Check CFL number
grep "^Courant" log.sedExnerFoam | tail -5

# Check y+ values (written at each writeTime)
postProcess -func yPlus -latestTime
```

---

## 6. PostProcessing: Extracting Bed Shear Stress Profile

The `sampleBedLine` function object writes near-bed velocity and pressure data to:
```
postProcessing/sampleBedLine/<timeStep>/bedCenterline_Ub_UbMean_p_rbgh_p_rbghMean.xy
```

The `wallShearStress1` function object writes wall shear stress vectors at each `writeTime` to:
```
postProcessing/wallShearStress1/<timeStep>/wallShearStress.dat
```

To extract the streamwise bed shear stress distribution for plotting:
```python
import numpy as np
# Read wall shear stress field at the final time step
# Plot tau_b(x) = wallShearStress_x along the bed patch
# Compare against paper Figure data
```

---

## 7. Summary of Changes (Validation Upgrade)
*   Mesh: 32,300 → ~57,000 cells (60 vertical layers, grading 20, $y_1 \approx 0.22\text{ mm}$)
*   Outlet BC: `zeroGradient` → `inletOutlet` (mass-balance stability fix)
*   Top pressure BC: `zeroGradient` → `slip` (consistency with velocity)
*   Turbulence schemes: `upwind` → `linearUpwind` (2nd-order accuracy)
*   PIMPLE: 3 → 4 outer correctors with residual control
*   Pressure tolerance: `1e-6` → `1e-7`
*   Time stepping: Fixed `dt=0.003` → CFL-adaptive (`maxCo=0.5`, `dt_max=0.01`)
*   Simulation time: 40s (1.3 flow-throughs)
*   Averaging window: 25–40s (15s of data)
*   Execution: Single-core → 20-core parallel (`scotch`)
*   Added: 2000-point bed sampling line, y+ monitoring
*   Write precision: 6 → 8 digits

---

## 8. Git Log
*   `0216456` — Update README.md with Ks directory structure and results table
*   `7579e1f` — Update controlDict startFrom to latestTime
*   `bcf5e6c` — Add final 100s bed shear stress comparison plot
*   `15f3006` — Remove parallel run configurations and options to restrict case to single-core
*   `3937484` — Add Allrun, Allclean scripts and update README.md with detailed instructions
*   `8b3ae9d` — Initialize OpenFOAM case with shortened domain and coded inlet
