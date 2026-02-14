/**
 * @file config_types.hpp
 * @brief Configuration types for the 1D Euler solver
 *
 * Strongly-typed configuration structures parsed from TOML input files.
 */

#ifndef EULER1D_CONFIG_CONFIG_TYPES_HPP
#define EULER1D_CONFIG_CONFIG_TYPES_HPP

#include "../core/types.hpp"
#include <string>
#include <vector>
#include <optional>

namespace euler1d {

// =============================================================================
// Enumerations for runtime-selectable options
// =============================================================================

/// Available numerical flux schemes
enum class FluxScheme {
    LLF,      ///< Local Lax-Friedrichs
    Rusanov,  ///< Rusanov (identical to LLF for scalar max wavespeed)
    HLL,      ///< Harten-Lax-van Leer
    HLLC,     ///< HLL with Contact restoration
    MoversLE  ///< MoversLE flux with adaptive dissipation
};

/// Available slope limiters for MUSCL reconstruction
enum class Limiter {
    None,      ///< No limiting (first order or unlimited)
    Minmod,    ///< Minmod limiter (most diffusive)
    VanLeer,   ///< Van Leer limiter
    Superbee,  ///< Superbee limiter (least diffusive)
    MC         ///< Monotonized Central limiter
};

/// Available time integration schemes
enum class TimeIntegrator {
    ExplicitEuler,  ///< Forward Euler (first order)
    SSPRK3          ///< Strong Stability Preserving RK3 (third order)
};

/// Available boundary condition types
enum class BoundaryType {
    Transmissive,  ///< Zero-gradient (outflow)
    Reflective,    ///< Solid wall (u = 0)
    Periodic       ///< Periodic boundaries
};

/// Available equation of state models
enum class EosModel {
    IdealGas  ///< Ideal gas with constant gamma
};

/// Available initial condition types
enum class InitialConditionType {
    PiecewiseConstant,        ///< Multiple constant regions
    ShockEntropyInteraction   ///< Shock + sinusoidal entropy wave
};

// =============================================================================
// Configuration Structures
// =============================================================================

/// Simulation metadata
struct SimulationConfig {
    std::string equations = "euler_1d";
    std::string test_name = "unnamed";
};

/// Mesh configuration
struct MeshConfig {
    Real xmin = 0.0;
    Real xmax = 1.0;
    int num_cells = 100;
};

/// Time stepping configuration
struct TimeConfig {
    Real cfl = 0.5;
    Real final_time = 1.0;
    TimeIntegrator integrator = TimeIntegrator::SSPRK3;
};

/// Numerical scheme configuration
struct NumericsConfig {
    int order = 1;  ///< 1 = first order, 2 = second order (MUSCL)
    FluxScheme flux = FluxScheme::LLF;
    Limiter limiter = Limiter::VanLeer;
};

/// Equation of state configuration
struct EosConfig {
    EosModel model = EosModel::IdealGas;
    Real gamma = 1.4;  ///< For ideal gas
};

/// Boundary conditions configuration
struct BoundaryConfig {
    BoundaryType left = BoundaryType::Transmissive;
    BoundaryType right = BoundaryType::Transmissive;
};

/// A constant region for piecewise initial conditions
struct Region {
    Real x_left = 0.0;
    Real x_right = 1.0;
    Real rho = 1.0;
    Real u = 0.0;
    Real p = 1.0;
};

/// State for shock-entropy interaction (constant part)
struct ConstantState {
    Real rho = 1.0;
    Real u = 0.0;
    Real p = 1.0;
};

/// Sinusoidal state for shock-entropy interaction
struct SinusoidalState {
    Real rho_base = 1.0;
    Real rho_amplitude = 0.0;
    Real rho_frequency = 0.0;
    bool use_pi = true;  ///< true: sin(freq*pi*x), false: sin(freq*x)
    Real u = 0.0;
    Real p = 1.0;
};

/// Initial condition configuration
struct InitialConditionConfig {
    InitialConditionType type = InitialConditionType::PiecewiseConstant;

    // For piecewise constant
    std::vector<Region> regions;

    // For shock-entropy interaction
    Real discontinuity_position = 0.0;
    ConstantState left_state;
    SinusoidalState right_state;
};

/// Complete configuration for the solver
struct Config {
    SimulationConfig simulation;
    MeshConfig mesh;
    TimeConfig time;
    NumericsConfig numerics;
    EosConfig eos;
    BoundaryConfig boundary;
    InitialConditionConfig initial_condition;
};

// =============================================================================
// String conversion utilities
// =============================================================================

/// Convert string to FluxScheme
FluxScheme parse_flux_scheme(const std::string& str);

/// Convert string to Limiter
Limiter parse_limiter(const std::string& str);

/// Convert string to TimeIntegrator
TimeIntegrator parse_time_integrator(const std::string& str);

/// Convert string to BoundaryType
BoundaryType parse_boundary_type(const std::string& str);

/// Convert string to EosModel
EosModel parse_eos_model(const std::string& str);

/// Convert string to InitialConditionType
InitialConditionType parse_initial_condition_type(const std::string& str);

}  // namespace euler1d

#endif  // EULER1D_CONFIG_CONFIG_TYPES_HPP
