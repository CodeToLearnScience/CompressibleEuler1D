# CompressibleEuler1D

A modular, high-performance 1D Compressible Euler solver using modern C++23.

## Features

- **No virtual dispatch**: Uses `std::variant` + `std::visit` for zero-overhead polymorphism
- **Runtime selectable**: Flux schemes, limiters, time integrators, boundary conditions
- **Compile-time selectable precision**: `float` or `double`
- **TOML configuration**: Human-readable input files
- **Comprehensive testing**: GoogleTest-based unit and integration tests
- **Multiple output formats**: CSV and VTK for visualization

## Building

### Requirements

- C++23 compatible compiler (GCC 13+, Clang 17+)
- CMake 3.25+

### Quick Start

```bash
# Configure (double precision, Release build)
cmake -B build -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build build -j

# Run tests
ctest --test-dir build --output-on-failure

# Run a simulation
./build/euler1d data/test_case1.toml
```

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `EULER1D_PRECISION` | `double` | Floating point precision (`double` or `float`) |
| `EULER1D_BUILD_TESTS` | `ON` | Build unit tests |

Example with float precision:
```bash
cmake -B build -DEULER1D_PRECISION=float -DCMAKE_BUILD_TYPE=Release
```

## Usage

```bash
./euler1d <config.toml> [output_dir]
```

## Configuration File Format

```toml
[simulation]
equations = "euler_1d"
test_name = "sod_shock_tube"

[mesh]
xmin = 0.0
xmax = 1.0
num_cells = 1000

[time]
cfl = 0.5
final_time = 0.2
time_integrator = "ssprk3"  # "euler" or "ssprk3"

[numerics]
order = 2          # 1 = first order, 2 = second order (MUSCL)
flux = "hllc"      # "llf", "rusanov", "hll", "hllc"
limiter = "vanleer" # "none", "minmod", "vanleer", "superbee", "mc"

[eos]
model = "ideal_gas"
gamma = 1.4

[boundary_conditions]
left = "transmissive"   # "transmissive", "reflective", "periodic"
right = "transmissive"

[initial_condition]
type = "piecewise_constant"

[[initial_condition.region]]
x_left = 0.0
x_right = 0.5
rho = 1.0
u = 0.0
p = 1.0

[[initial_condition.region]]
x_left = 0.5
x_right = 1.0
rho = 0.125
u = 0.0
p = 0.1
```

## Numerical Schemes

### Flux Schemes

| Scheme | Description |
|--------|-------------|
| LLF | Local Lax-Friedrichs (most diffusive) |
| Rusanov | Rusanov flux (same as LLF) |
| HLL | Harten-Lax-van Leer (two-wave solver) |
| HLLC | HLL with Contact restoration (recommended) |

### Limiters

| Limiter | Description |
|---------|-------------|
| none | No limiting (first order) |
| minmod | Most diffusive TVD limiter |
| vanleer | Smooth TVD limiter (recommended) |
| superbee | Least diffusive TVD limiter |
| mc | Monotonized Central |

### Time Integrators

| Integrator | Order | Description |
|------------|-------|-------------|
| euler | 1 | Forward Euler (not recommended for production) |
| ssprk3 | 3 | Strong Stability Preserving RK3 (recommended) |

## Extending the Solver

### Adding a New Flux Scheme

1. Define struct in `include/euler1d/flux/flux.hpp`:
```cpp
struct MyFlux {
    template <typename Eos>
    ConservativeVars operator()(const ConservativeVars& U_L,
                                 const ConservativeVars& U_R,
                                 const Eos& eos) const noexcept {
        // Implementation
    }
};
```

2. Add to `FluxVariant`:
```cpp
using FluxVariant = std::variant<LLFFlux, RusanovFlux, HLLFlux, HLLCFlux, MyFlux>;
```

3. Add enum value and parser support in `config_types.hpp` and `parser.cpp`

4. Add factory case in `factory.cpp`

### Adding a New Limiter

Same pattern as flux schemes - define struct with `operator()(Real r)`, add to variant, update parser and factory.

## Output Files

- **CSV**: `test_name.csv` - columns: x, rho, u, p, E
- **VTK**: `test_name.vtk` - ParaView compatible structured grid

## Test Cases

The `data/` directory contains 12 standard test cases:

| # | Name | Key Feature |
|---|------|-------------|
| 1 | Toro shock tube | Sonic point |
| 2 | Toro overheating | Colliding flows |
| 3 | Toro strong shock | Pressure ratio 100000:1 |
| 4 | Toro several shocks | High velocities |
| 5 | Toro slowly moving contact | Moving frame |
| 6 | Steady contact | Density jump only |
| 7 | Zhang & Shu Mach 2 | Long time evolution |
| 8 | Slowly moving shock | Weak shock |
| 9 | Slowly moving contact | Advected contact |
| 10 | Two discontinuities | Blast wave interaction |
| 11 | Balsara & Shu | Shock-entropy (5πx) |
| 12 | Shu & Osher | Shock-entropy (5x) |

## License

See LICENSE file.
