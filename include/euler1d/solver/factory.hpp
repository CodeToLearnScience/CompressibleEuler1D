/**
 * @file factory.hpp
 * @brief Factory functions for creating solver components
 */

#ifndef EULER1D_SOLVER_FACTORY_HPP
#define EULER1D_SOLVER_FACTORY_HPP

#include "../config/config_types.hpp"
#include "../eos/eos.hpp"
#include "../flux/flux.hpp"
#include "../reconstruction/limiter.hpp"
#include "../boundary/boundary.hpp"
#include "../time/time_integrator.hpp"
#include "../initial/initial_condition.hpp"

namespace euler1d {

/// Create EOS from config
EosVariant create_eos(const EosConfig& config);

/// Create flux scheme from enum
FluxVariant create_flux(FluxScheme scheme);

/// Create limiter from enum
LimiterVariant create_limiter(Limiter lim);

/// Create boundary condition from enum
BoundaryVariant create_boundary(BoundaryType type);

/// Create time integrator from enum
TimeIntegratorVariant create_time_integrator(TimeIntegrator integ);

}  // namespace euler1d

#endif  // EULER1D_SOLVER_FACTORY_HPP
