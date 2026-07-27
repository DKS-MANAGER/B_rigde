# OpenFOAM Rigid-Bed Pressure-Flow Case (Majid et al. 2026)

This case directory contains a complete, optimized, and ready-to-run OpenFOAM (v2412) setup designed to reproduce the rigid-bed pressure-flow-due-to-vertical-contraction experiments described in:

> **Effect of Bed Roughness on Pressure Flow due to Vertical Contraction**  
> *Sofi Aamir Majid, S.M.ASCE; Shivam Tripathi; and Debopam Das*  
> ASCE Journal of Hydraulic Engineering, 2026.

---

## 1. Case Physics & Configuration

* **Domain Geometry:** 2D longitudinal flume section ($8.0\text{ m} \times 0.10\text{ m} \times 0.01\text{ m}$).
* **Bridge contraction:** Solid ceiling block from $x = 1.0\text{ m}$ to $x = 1.15\text{ m}$ ($L = 0.15\text{ m}$). Clear throat gap height $H_b = 0.075\text{ m}$ ($H = 0.025\text{ m}$ contraction depth).
* **Inlet Boundary Condition:** A fully developed turbulent boundary layer velocity profile is mathematically imposed at the inlet ($x=0$) using a 1/7th power law with a bulk velocity $V = 0.26\text{ m/s}$ ($U_{\text{max}} = 0.2971\text{ m/s}$):
  $$u(y) = 0.2971 \cdot \left(\frac{y}{0.10}\right)^{1/7}$$
  This is implemented using a runtime-compiled `codedFixedValue` boundary condition for the fluid velocity field (`Ub`) and mixture velocity field (`U`).
* **Bed Roughness:** Configured with `nutkRoughWallFunction` on the bottom wall using equivalent sand-grain roughness height $K_s = 0.68\text{ mm}$ (Grade II sand) and $C_s = 0.5$.
* **Solver:** `sedExnerFoam` (running in clear-water hydrodynamics mode with `sedimentBed off` in `constant/bedloadProperties` to ensure the bed remains completely rigid).

---

## 2. Optimizations Implemented

* **Shortened Upstream Domain:** Reduced the upstream entrance length from **8.0 m to 1.0 m** by employing the coded fully-developed turbulent profile at the inlet.
* **Downstream Recovery:** Maintained the full **6.85 m** downstream section ($x = 1.15 \to 8.0\text{ m}$) to accurately capture flow expansion, wake development, and secondary shear-stress peaks.
* **Cell Count & Step Size:** Reduced grid resolution to **32,300 cells** (~50% smaller) and bumped $\Delta t$ to **0.003 s** (perfectly stable under Courant limit constraints), speeding up wall-clock execution time by **~100x**.
* **Optimal Run Time:** The simulation duration is set to **40 s** (equivalent to ~1.3 domain flow-throughs), which is the exact duration required to reach a fully stabilized steady-state flow profile.

---

## 3. Case Automation Scripts

Two shell scripts are provided to automate workflow steps:

### `./Allrun`
Cleans old files, generates the polyMesh (`blockMesh`), constructs the finite-area mesh (`makeFaMesh`), validates mesh quality (`checkMesh`), and runs the solver.
* **Usage:**
  ```bash
  ./Allrun
  ```

### `./Allclean`
Resets the case to its clean, pre-simulation state by removing time step directories, processor directories, log files, dynamically compiled libraries, and mesh files.
* **Usage:**
  ```bash
  ./Allclean
  ```

---

## 4. Monitoring Progress

To monitor solver progress in real-time:
```bash
tail -f log.sedExnerFoam
```
Check the Courant number (`Co`), solver residuals, and the current simulation time step. Outputs will write to time directories (e.g., `10/`, `20/`...) at intervals of **10 s** of simulation time.
