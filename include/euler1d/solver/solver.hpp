/**
 * @file solver.hpp
 * @brief Main solver orchestration for the 1D Euler equations
 */

#ifndef EULER1D_SOLVER_SOLVER_HPP
#define EULER1D_SOLVER_SOLVER_HPP

#include "../core/types.hpp"
#include "../config/config_types.hpp"
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"
#include "../flux/flux.hpp"
#include "../reconstruction/limiter.hpp"
#include "../reconstruction/muscl.hpp"
#include "../boundary/boundary.hpp"
#include "../initial/initial_condition.hpp"
#include "../time/time_integrator.hpp"
#include <functional>
#include <iostream>
#include <memory>

namespace euler1d {

/**
 * @brief Main solver class for 1D Euler equations
 *
 * Orchestrates mesh, EOS, flux scheme, reconstruction, boundaries,
 * time integration, and solution output.
 */
class Solver {
public:
    /// Construct solver from configuration
    explicit Solver(const Config& config);

    /// Run simulation to final time
    void run();

    /// Get current solution (conservative variables)
    [[nodiscard]] const ConservativeArray& solution() const noexcept { return U_; }

    /// Get mesh
    [[nodiscard]] const Mesh1D& mesh() const noexcept { return mesh_; }

    /// Get current simulation time
    [[nodiscard]] Real time() const noexcept { return time_; }

    /// Get test name from config
    [[nodiscard]] const std::string& test_name() const noexcept { return config_.simulation.test_name; }

    /// Convert solution to primitive variables
    [[nodiscard]] PrimitiveArray to_primitive() const;

private:
    /// Compute RHS: dU/dt = -d(F)/dx
    void compute_rhs(std::span<const ConservativeVars> U, std::span<ConservativeVars> dU);

    /// Compute stable timestep based on CFL condition
    [[nodiscard]] Real compute_dt() const;

    /// Apply boundary conditions
    void apply_boundaries();

    /// Convert solution to primitive variables (internal buffer)
    void update_primitives();

    Config config_;
    Mesh1D mesh_;
    EosVariant eos_;
    FluxVariant flux_;
    LimiterVariant limiter_;
    BoundaryVariant bc_left_;
    BoundaryVariant bc_right_;
    TimeIntegratorVariant time_integrator_;
    InitialConditionVariant initial_condition_;

    ConservativeArray U_;       ///< Current solution (conservative)
    PrimitiveArray W_;          ///< Current solution (primitive)
    ConservativeArray fluxes_;  ///< Interface fluxes

    Real time_ = 0;
    int order_ = 1;
};

}  // namespace euler1d

#endif  // EULER1D_SOLVER_SOLVER_HPP
