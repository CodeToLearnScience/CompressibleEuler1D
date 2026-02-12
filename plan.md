# Plan: 1D Compressible Euler Finite Volume Solver

A modular, high-performance 1D Euler solver using modern C++23 with compile-time polymorphism and runtime configuration via TOML. The design prioritizes extensibility (new fluxes, limiters, time integrators) without virtual dispatch.

## High-Level Architecture

```
CompressibleEuler1D/
├── CMakeLists.txt
├── cmake/
│   └── dependencies.cmake
├── include/euler1d/
│   ├── config/          # TOML parsing & configuration types
│   ├── core/            # Precision typedef, constants, conservative/primitive vars
│   ├── eos/             # Equation of state (ideal gas, extensible)
│   ├── mesh/            # 1D mesh with ghost cells
│   ├── flux/            # LLF, Rusanov, HLL, HLLC
│   ├── reconstruction/  # MUSCL, limiters
│   ├── boundary/        # Boundary conditions
│   ├── initial/         # Initial condition generators
│   ├── time/            # Euler, SSPRK3
│   ├── solver/          # Main solver orchestration
│   └── io/              # Output (CSV, VTK)
├── src/                 # Implementation files
├── tests/               # GTest unit & integration tests
├── apps/
│   └── main.cpp         # Entry point
└── docs/
    └── README.md
```

---

## Steps

### 1. Project Setup & Build System

1. Create root `CMakeLists.txt` with:
   - C++23 standard requirement (`CMAKE_CXX_STANDARD 23`)
   - Compile-time precision option: `-DPRECISION=double|float`
   - FetchContent for dependencies: **toml++**, **GoogleTest**
   - Targets: `euler1d_lib` (static library), `euler1d` (executable), `euler1d_tests`

2. Create `cmake/dependencies.cmake` for FetchContent declarations

### 2. Core Types & Precision

3. Create `include/euler1d/core/types.hpp`:
   - `using Real = PRECISION;` (compile-time selectable)
   - `struct ConservativeVars { Real rho, rho_u, E; }`
   - `struct PrimitiveVars { Real rho, u, p; }`
   - Conversion functions: `conservative_to_primitive()`, `primitive_to_conservative()`
   - Template on `Real` throughout

4. Create `include/euler1d/core/constants.hpp`:
   - Physical and numerical constants

### 3. Equation of State (EOS)

5. Create `include/euler1d/eos/eos.hpp`:
   - `struct IdealGas { Real gamma; }` with methods: `pressure()`, `sound_speed()`, `internal_energy()`, `enthalpy()`
   - Use `std::variant<IdealGas, ...>` for extensibility to other EOS models
   - Free functions dispatched via `std::visit`

### 4. TOML Configuration Parser

6. Create `include/euler1d/config/config_types.hpp`:
   - Strongly-typed structs for each TOML section:
     - `SimulationConfig`, `MeshConfig`, `TimeConfig`, `NumericsConfig`, `EosConfig`, `BoundaryConfig`, `InitialConditionConfig`
   - Enums: `FluxScheme`, `Limiter`, `TimeIntegrator`, `BoundaryType`, `InitialConditionType`
   - `[[initial_condition.region]]` as `std::vector<Region>`
   - Special handling for `shock_entropy_interaction` with sinusoidal parameters

7. Create `include/euler1d/config/parser.hpp` and `src/config/parser.cpp`:
   - `Config parse_config(const std::filesystem::path& file)`
   - Validate required fields, throw descriptive errors
   - Support all 12 test case formats

### 5. 1D Mesh

8. Create `include/euler1d/mesh/mesh.hpp`:
   - `class Mesh1D`: stores `xmin`, `xmax`, `num_cells`, `dx`
   - Ghost cell support (2 ghost cells each side for 2nd order)
   - Cell center coordinates: `x(i)` accessor
   - Index helpers: `first_interior()`, `last_interior()`, `num_ghosts()`

### 6. Flux Schemes

9. Create `include/euler1d/flux/flux.hpp`:
   - `using FluxVariant = std::variant<LLF, Rusanov, HLL, HLLC>;`
   - Each flux as a struct with `operator()(left, right, eos) -> ConservativeVars`

10. Implement flux schemes in `src/flux/`:
    - `llf.cpp`: Local Lax-Friedrichs
    - `rusanov.cpp`: Rusanov (similar to LLF)
    - `hll.cpp`: HLL two-wave solver
    - `hllc.cpp`: HLLC with contact restoration

### 7. Reconstruction & Limiters

11. Create `include/euler1d/reconstruction/limiter.hpp`:
    - `using LimiterVariant = std::variant<NoLimiter, Minmod, VanLeer, Superbee, MC>;`
    - Each limiter: `operator()(r) -> Real` (slope limiter function)

12. Create `include/euler1d/reconstruction/muscl.hpp`:
    - `template<typename Limiter> reconstruct(U, i, limiter) -> (U_L, U_R)`
    - Characteristic-based or component-wise reconstruction option
    - First order (no reconstruction) and second order (MUSCL)

### 8. Boundary Conditions

13. Create `include/euler1d/boundary/boundary.hpp`:
    - `using BoundaryVariant = std::variant<Transmissive, Reflective, Periodic, ...>;`
    - `apply_left(U_array, mesh)`, `apply_right(U_array, mesh)`
    - Transmissive: copy interior to ghost cells

### 9. Initial Conditions

14. Create `include/euler1d/initial/initial_condition.hpp`:
    - `using ICVariant = std::variant<PiecewiseConstant, ShockEntropyInteraction>;`
    - `PiecewiseConstant`: loop over regions, set values based on `x_left`, `x_right`
    - `ShockEntropyInteraction`: 
      - Left of discontinuity: constant state
      - Right: `rho = rho_base + rho_amplitude * sin(rho_frequency * [pi*x or x])`
      - Parse `rho_function = "pi"` vs `"linear"`

### 10. Time Integration

15. Create `include/euler1d/time/time_integrator.hpp`:
    - `using TimeIntegratorVariant = std::variant<ExplicitEuler, SSPRK3>;`
    - Each integrator: `advance(U, dt, rhs_function) -> U_new`

16. Implement in `src/time/`:
    - `explicit_euler.cpp`: `U^{n+1} = U^n + dt * L(U^n)`
    - `ssprk3.cpp`: 3-stage SSP-RK3:
      ```
      U^(1) = U^n + dt*L(U^n)
      U^(2) = 3/4*U^n + 1/4*(U^(1) + dt*L(U^(1)))
      U^(3) = 1/3*U^n + 2/3*(U^(2) + dt*L(U^(2)))
      ```

### 11. Solver Orchestration

17. Create `include/euler1d/solver/solver.hpp` and `src/solver/solver.cpp`:
    - `class Solver` holding all variant types:
      ```cpp
      FluxVariant flux_;
      LimiterVariant limiter_;
      TimeIntegratorVariant time_integrator_;
      BoundaryVariant bc_left_, bc_right_;
      EosVariant eos_;
      ```
    - `compute_rhs()`: reconstruction → flux computation → flux differencing
    - `compute_dt()`: CFL-based timestep
    - `run()`: main time loop with progress output

18. Create factory functions in `src/solver/factory.cpp`:
    - `create_flux(FluxScheme) -> FluxVariant`
    - `create_limiter(Limiter) -> LimiterVariant`
    - etc.

### 12. Output (I/O)

19. Create `include/euler1d/io/output.hpp` and `src/io/`:
    - `csv_writer.cpp`: `write_csv(path, mesh, U, time)` → columns: x, rho, u, p, E
    - `vtk_writer.cpp`: `write_vtk(path, mesh, U, time)` → VTK legacy format for ParaView

### 13. Main Application

20. Create `apps/main.cpp`:
    - Parse command-line argument for TOML file path
    - Load config → create solver → run → write output
    - Print summary (time, cells, final time, etc.)

### 14. Testing with GTest

21. Create test structure under `tests/`:
    - `tests/CMakeLists.txt`
    - `tests/test_eos.cpp`: ideal gas thermodynamic relations
    - `tests/test_config_parser.cpp`: parse all 12 test cases
    - `tests/test_flux.cpp`: flux consistency (symmetry, conservation)
    - `tests/test_limiter.cpp`: limiter bounds and TVD properties
    - `tests/test_reconstruction.cpp`: accuracy on smooth data
    - `tests/test_boundary.cpp`: ghost cell values
    - `tests/test_initial_condition.cpp`: piecewise & sinusoidal
    - `tests/test_time_integrator.cpp`: order of accuracy
    - `tests/test_solver_integration.cpp`: Sod shock tube convergence

### 15. Documentation

22. Create `docs/README.md`:
    - Build instructions, usage, TOML format reference
    - Algorithm descriptions
    - Extending the solver (adding new flux/limiter)

23. Add Doxygen-style comments to all public headers

---

## Verification

1. **Build test**: `cmake -B build -DCMAKE_BUILD_TYPE=Release -DPRECISION=double && cmake --build build`
2. **Unit tests**: `ctest --test-dir build --output-on-failure`
3. **Integration**: Run all 12 test cases:
   ```bash
   for f in data/*.toml; do ./build/apps/euler1d "$f"; done
   ```
4. **Validation**: Compare Sod shock tube (test_case1) against analytical solution
5. **Performance**: Profile with `perf` or timing output

---

## Decisions

- **std::variant + std::visit** over virtual dispatch: zero-overhead polymorphism, type-safe, compiler can inline and optimize
- **C++23**: enables `std::print`, improved ranges, `[[assume]]` for optimization hints
- **MUSCL reconstruction on primitive variables**: more robust than conservatives for Euler
- **2 ghost cells**: supports up to 4th-order stencils for future extension
- **Component-wise limiting first**: simpler, characteristic-based as optional enhancement
- **Factory pattern with enums**: clean mapping from config strings to variant types
