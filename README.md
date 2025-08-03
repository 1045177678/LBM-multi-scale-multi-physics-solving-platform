# LBM-multi-scale-multi-physics-solving-platform
A general multiscale multiphysics solver based on the lattice Boltzmann method
# LBM_newSolver — 2-D Multi-Level LBM (D2Q9, BGK/MRT, Scalar Fields, r² Subcycling)

**LBM_newSolver** is a small, extensible C++17 framework for 2-D multi-level Lattice Boltzmann (LBM) using D2Q9. It supports flow and scalar fields, time-varying operators, SI↔LU unit conversion per level, and r² time subcycling between coarse/fine grids with 6-moment restriction/prolongation. The code builds directly with Visual Studio 2022; no CMake required (optional).

> For peer review: this repo includes a minimal runnable example (driven cavity), config schema, and troubleshooting notes.

---

## Features

- **Multi-level grids + r² subcycling**  
  If `dx_c/dx_f = r` (integer), then `dt_f = dt_c / r²`. Fine levels advance `r²` steps per coarse step with synchronization points.

- **Unified fields & methods**  
  Flow (`kind: "flow"`) and scalar (`kind: "scalar"`) both use D2Q9 with BGK/MRT backends. Scalars can be advected by a chosen flow field.

- **6-moment transfer between levels**  
  On overlaps, D2Q9 populations are projected to 6 hydrodynamic moments for restriction/prolongation to avoid spurious oscillations; then reconstructed.

- **Time-varying operators (TimeExpr)**  
  Boundary conditions, body forces, heat sources etc. accept constant / table / small expression forms (`sin`, `cos`, `exp`, `log`, `min`, `max`, `pi`, `t`).

- **Layerwise SI↔LU unit system**  
  Each level declares `CL` (m per cell) and optionally `T0` (s per step). Given `tau` & `nu_phys`, solver back-computes consistent `T0` per level so that all levels keep an optimal `tau` while time steps remain integer-related.

- **Adaptive output**  
  Tecplot FEQUADRILATERAL/ BLOCK `.dat` with only the variables of active fields (e.g., `rho, ux, uy, phi`), plus timestamps.
  - **Boundary conditions as Operators**  
  In this framework, all boundary conditions (e.g. velocity or no-slip BCs) are implemented as `Operator` instances. During each time step, the engine automatically loops over every operator targeting the current level and calls `op->apply(ctx)`. This unified “operators” stage ensures that boundary conditions are applied seamlessly alongside body-forces, sources, and other custom operators before and/or after the collision-streaming phases, without any special casing in the core time‐step logic.

---

## Directory Layout (suggested)

```
.
├─ src/
│  ├─ main.cpp
│  ├─ lbm_core.h
│  ├─ lbm_core.cpp
│  ├─ lbm_tecout.cpp
│  └─ json.hpp
├─ examples/
│  ├─ mesh_meta_cell.json
│  └─ config.json
├─ readme       # this file
└─ LICENSE
```

---

## Build & Run

### Visual Studio 2022 (Windows)
1. Create an empty x64 C++ project and add all files in `src/`.  
2. Set **C++ Language Standard** to C++17.  
3. Provide **nlohmann/json** either as a single header (`json.hpp`) or via NuGet/vcpkg.  
4. In **Project → Debugging**, set:
   - **Command Arguments**: `examples/config.json`  
   - **Working Directory**: `$(ProjectDir)`  
5. Build & run.

### Command Line
```bash
g++ src/*.cpp -std=c++17 -O3 -o lbm
./lbm examples/config.json
```

---

## Minimal Example Files

- **mesh_meta_cell.json**: grid metadata (levels, dims, spacing).  
- **config.json**: solver steps, units, fields, operators.

---

## License

MIT. Please cite this repository when publishing results.

