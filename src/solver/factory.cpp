/**
 * @file factory.cpp
 * @brief Factory functions implementation
 */

#include "euler1d/solver/factory.hpp"

namespace euler1d {

/// Create EOS from config
EosVariant create_eos(const EosConfig& config) {
    switch (config.model) {
        case EosModel::IdealGas:
            return IdealGas{config.gamma};
    }
    return IdealGas{};
}

/// Create flux scheme from config
FluxVariant create_flux(FluxScheme scheme) {
    switch (scheme) {
        case FluxScheme::LLF: return LLFFlux{};
        case FluxScheme::Rusanov: return RusanovFlux{};
        case FluxScheme::HLL: return HLLFlux{};
        case FluxScheme::HLLC: return HLLCFlux{};
    }
    return LLFFlux{};
}

/// Create limiter from config
LimiterVariant create_limiter(Limiter lim) {
    switch (lim) {
        case Limiter::None: return NoLimiter{};
        case Limiter::Minmod: return MinmodLimiter{};
        case Limiter::VanLeer: return VanLeerLimiter{};
        case Limiter::Superbee: return SuperbeeLimiter{};
        case Limiter::MC: return MCLimiter{};
    }
    return VanLeerLimiter{};
}

/// Create boundary condition from config
BoundaryVariant create_boundary(BoundaryType type) {
    switch (type) {
        case BoundaryType::Transmissive: return TransmissiveBoundary{};
        case BoundaryType::Reflective: return ReflectiveBoundary{};
        case BoundaryType::Periodic: return PeriodicBoundary{};
    }
    return TransmissiveBoundary{};
}

/// Create time integrator from config
TimeIntegratorVariant create_time_integrator(TimeIntegrator integ) {
    switch (integ) {
        case TimeIntegrator::ExplicitEuler: return ExplicitEuler{};
        case TimeIntegrator::SSPRK3: return SSPRK3{};
    }
    return SSPRK3{};
}

}  // namespace euler1d
