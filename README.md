# CFD Validation Report: Pressure-Flow due to Vertical Bridge Contraction

This repository contains the OpenFOAM (v2412) implementation, grid setup, validation profiles, and post-processing tools designed to reproduce the rigid-bed pressure-flow-due-to-vertical-contraction flume experiments reported by:

> **Effect of Bed Roughness on Pressure Flow due to Vertical Contraction**  
> *Sofi Aamir Majid, S.M.ASCE; Shivam Tripathi; and Debopam Das*  
> **Journal of Hydraulic Engineering, ASCE (Volume 152, Issue 3, January 2026)**  
> **DOI:** [10.1061/JHEND8.HYENG-14490](https://doi.org/10.1061/JHEND8.HYENG-14490)

---

## 1. Project Title & Overview
This project validates the bed shear stress ($\tau_b$) distribution along a flume channel bed subjected to a sudden vertical constriction. Using a customized OpenFOAM hydrodynamic solver (`sedExnerFoam` running in hydrodynamics-only mode), we evaluate three different bed roughness grades ($K_s = 0.33\text{ mm}$, $0.68\text{ mm}$, and $1.90\text{ mm}$). 

The model resolves the wall shear stress signature, validating the numerical output directly against the PIV experimental measurements from the ASCE 2026 publication.

---

## 2. Hydraulic & Physical Parameters
The numerical setup models a 2D slice of the experimental hydraulic flume. The exact physical and flow properties are detailed below:

| Parameter | Symbol | Value | Unit | Description |
| :--- | :---: | :---: | :---: | :--- |
| **Approach Flow Depth** | $H_a$ | $0.10$ | $\text{m}$ | Constant water level upstream |
| **Channel Width** | $Z_w$ | $0.30$ | $\text{m}$ | Transverse width of the flume |
| **Channel Bed Slope** | $S_0$ | $0.018$ | $\%$ | Streamwise slope ($0.00018\text{ m/m}$) |
| **Fluid Density** | $\rho_f$ | $1000$ | $\text{kg/m}^3$ | Density of water at $20^\circ\text{C}$ |
| **Kinematic Viscosity** | $\nu_f$ | $1.0 \times 10^{-6}$ | $\text{m}^2/\text{s}$ | Viscosity of water |
| **Bulk Flow Velocity** | $V_a$ | $0.26$ | $\text{m/s}$ | Average approach flow rate velocity |
| **Grade I Roughness** | $K_{s,1}$ | $0.33$ | $\text{mm}$ | Median diameter $d_{50}$ for Grade I sand |
| **Grade II Roughness** | $K_{s,2}$ | $0.68$ | $\text{mm}$ | Median diameter $d_{50}$ for Grade II sand |
| **Grade III Roughness** | $K_{s,3}$ | $1.90$ | $\text{mm}$ | Median diameter $d_{50}$ for Grade III sand |

---

## 3. Repository & Directory Map
The workspace is structured into three self-contained OpenFOAM case folders representing each bed roughness grade:

```directory
E:\DKS\B_ridgi
├── Ks_0.33/                        # Case directory for Grade I (Ks = 0.33 mm)
│   ├── 0/                          # Boundary and initial conditions (t=0)
│   ├── constant/                   # Transport & turbulence properties, faMesh
│   ├── system/                     # blockMeshDict, fvSchemes, fvSolution
│   └── Allrun, Allclean            # Automation run/clean scripts
├── Ks_1.9/                         # Case directory for Grade III (Ks = 1.90 mm)
│   ├── 0/, constant/, system/      # Case settings matching Grade III
│   └── Allrun, Allclean            # Automation run/clean scripts
├── 0/                              # Case settings for Grade II (Ks = 0.68 mm) - Root Case
├── constant/                       # Transport & mesh properties for Root Case
├── system/                         # blockMeshDict, solvers, controlDict for Root Case
├── Allrun, Allclean                # Automation run/clean scripts for Root Case
├── plot_comparison.py              # Python utility to parse time steps and plot validation curve
└── bed_shear_stress_comparison.png # Generated validation comparison plot
```

### Key Dictionary Configurations:
*   **`system/blockMeshDict`**: Generates a structured multi-block mesh spanning $x \in [0.0\text{ m}, 8.0\text{ m}]$ and $y \in [0.0\text{ m}, 0.1\text{ m}]$. The vertical bridge contraction is located at $x \in [1.0\text{ m}, 1.15\text{ m}]$ with a ceiling height of $y = 0.075\text{ m}$ (representing a 25% height constriction).
*   **`constant/transportProperties`**: Defines fluid properties (density, kinematic viscosity).
*   **`system/fvSchemes`**: Controls spatial and temporal discretization. Modified to use stable bounded schemes.
*   **`system/fvSolution`**: Dictates linear solvers (PCG, PBiCGStab), tolerances, and PIMPLE outer loop loops.

---

## 4. Solver Fixes & Stability Tuning Log

This validation campaign resolved two major numerical issues that initially affected the simulation accuracy:

### 4.1 Fuhrman Wall Function Parameter Correction (Overlapping Curves Fix)
*   **The Issue:** Previously, the bed shear stress profiles for all three roughness grades ($0.33\text{ mm}$, $0.68\text{ mm}$, and $1.90\text{ mm}$) overlapped, outputting identical results.
*   **The Discovery:** The custom compiled library `libroughWallFunctions.so` uses the Menter-Esch/Fuhrman formulation (`fuhrmanOmegaWallFunction`) for $\omega$. The C++ source code expects the equivalent roughness parameter named **`kn`** (Nikuradse roughness). Because the dictionaries initially provided the standard parameter **`Ks`**, it was ignored, and the wall function silently fell back to its default value of `kn = 1e-6` in all three cases.
*   **The Fix:** Updated the boundary conditions on the `bed` patch in `0/omega` and `0/omega.b` to supply the correct parameter:
    ```cpp
    bed
    {
        type            fuhrmanOmegaWallFunction;
        kn              0.00033; // 0.00068 for Grade II, 0.0019 for Grade III
        value           uniform 4.153;
    }
    ```
    This successfully activated the roughness boundary shift, leading to the physical separation of the shear stress profiles.

### 4.2 Discretization Scheme Correction (Fluctuations Fix)
*   **The Issue:** The resolved shear stress profiles showed high-frequency spatial fluctuations (wiggles) along the bed centerline, especially downstream of the contraction block.
*   **The Discovery:** The convection terms for the turbulent parameters (`k` and `omega`) were using the second-order `linearUpwind` scheme, which triggered localized numerical oscillations near the high-gradient wall boundaries on this fine structured grid.
*   **The Fix:** Changed the convection schemes of `k` and `omega` in `system/fvSchemes` to the first-order bounded **`Gauss upwind`** scheme:
    ```cpp
    div(phi,k)      Gauss upwind;
    div(phi,omega)  Gauss upwind;
    ```
    This completely eliminated the spatial oscillations, resulting in smooth, physically realistic curves.

---

## 5. Prerequisites & Environment Setup
To run the simulations, the following environment is required:
*   **OpenFOAM v2412** (installed and sourced in your terminal/WSL).
*   **custom compiled libraries:** `libroughWallFunctions.so` must be present in the user library directory (`$FOAM_USER_LIBBIN`).
*   **Python 3** with `matplotlib` and `numpy` installed for validation plotting.

Source the OpenFOAM environment in your shell before running:
```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
```

---

## 6. Step-by-Step Execution Guide

To run a simulation case (e.g., Grade I):
```bash
# Navigate to case directory
cd Ks_0.33

# Run the automated execution script
./Allrun
```

The `./Allrun` script automates the following sequence:
1.  **Cleanup:** Cleans previous run logs and time steps (`Allclean`).
2.  **Mesh Generation:** Runs `blockMesh` to generate the 57,000-cell block grid.
3.  **Finite-Area Mesh:** Runs `makeFaMesh` for boundary zone definitions.
4.  **Parallel Decomposition:** Decomposes the mesh into 4 processors (`decomposePar`).
5.  **Parallel Launch:** Executes the solver in parallel:
    ```bash
    mpirun -np 4 sedExnerFoam -parallel > log.sedExnerFoam 2>&1 &
    ```

To reconstruct the final parallel time step folders for post-processing/ParaView:
```bash
reconstructPar -time 40
```

To reconstruct all written time steps:
```bash
reconstructPar
```

---

## 7. Post-Processing & Validation Metrics

The streamwise distribution of the bed shear stress ($\tau_b$) is extracted from the `wallShearStress` vector field:
$$\tau_b = \rho_f \times |\tau_{w,x}|$$
*(where OpenFOAM outputs kinematic shear stress $\tau_{w} / \rho_f$, requiring multiplication by density $\rho_f = 1000\text{ kg/m}^3$ to obtain Pascal values).*

A Python script `plot_comparison.py` is included in the root folder. Running this script automatically parses the final time step folder ($t = 40\text{ s}$), extracts the bed centerline profile, and plots the CFD results against the ASCE 2026 experimental points:
```bash
python3 plot_comparison.py
```
This updates the validation figure **`bed_shear_stress_comparison.png`** in the root directory.

---

## 8. Acknowledgements & References
*   **Original Publication:** Majid et al., "Effect of Bed Roughness on Pressure Flow due to Vertical Contraction", *Journal of Hydraulic Engineering, ASCE*, 2026. DOI: [10.1061/JHEND8.HYENG-14490](https://doi.org/10.1061/JHEND8.HYENG-14490).
*   **Solver and Library Support:** Development and integration of the custom rough wall functions library is credited to the **LEGI** (Laboratoire des Ecoulements Geophysiques et Industriels) SedFoam research community.
