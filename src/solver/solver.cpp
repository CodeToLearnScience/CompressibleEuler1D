/**
 * @file solver.cpp
 * @brief Solver implementation
 */

#include "euler1d/solver/solver.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <format>
#include <print>

namespace euler1d {

namespace {

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

}  // namespace

Solver::Solver(const Config& config)
    : config_{config},
      mesh_{config.mesh.xmin, config.mesh.xmax, config.mesh.num_cells},
      eos_{create_eos(config.eos)},
      flux_{create_flux(config.numerics.flux)},
      limiter_{create_limiter(config.numerics.limiter)},
      bc_left_{create_boundary(config.boundary.left)},
      bc_right_{create_boundary(config.boundary.right)},
      time_integrator_{create_time_integrator(config.time.integrator)},
      initial_condition_{create_initial_condition(config.initial_condition)},
      order_{config.numerics.order} {

    // Allocate solution arrays
    const auto n = static_cast<std::size_t>(mesh_.total_cells());
    U_.resize(n);
    W_.resize(n);
    fluxes_.resize(n + 1);  // n+1 interfaces

    // Apply initial condition
    std::visit([this](const auto& eos) {
        std::visit([this, &eos](const auto& ic) {
            ic.apply(U_, mesh_, eos);
        }, initial_condition_);
    }, eos_);

    // Apply boundary conditions
    apply_boundaries();

    // Initialize primitives
    update_primitives();
}

void Solver::apply_boundaries() {
    apply_left_boundary(bc_left_, U_, mesh_);
    apply_right_boundary(bc_right_, U_, mesh_);
}

void Solver::update_primitives() {
    std::visit([this](const auto& eos) {
        for (std::size_t i = 0; i < U_.size(); ++i) {
            W_[i] = eos.to_primitive(U_[i]);
        }
    }, eos_);
}

Real Solver::compute_dt() const {
    Real max_speed = Real{0};

    std::visit([this, &max_speed](const auto& eos) {
        for (int i = mesh_.first_interior(); i <= mesh_.last_interior(); ++i) {
            const auto& U = U_[static_cast<std::size_t>(i)];
            const Real u = std::abs(U.rho_u / U.rho);
            const Real c = eos.sound_speed(U);
            max_speed = std::max(max_speed, u + c);
        }
    }, eos_);

    if (max_speed < constants::epsilon) {
        max_speed = Real{1};  // Avoid division by zero
    }

    return config_.time.cfl * mesh_.dx() / max_speed;
}

void Solver::compute_rhs(std::span<const ConservativeVars> U, std::span<ConservativeVars> dU) {
    // Update primitives from U
    std::visit([&U, this](const auto& eos) {
        for (std::size_t i = 0; i < U.size(); ++i) {
            W_[i] = eos.to_primitive(U[i]);
        }
    }, eos_);

    // Compute fluxes at each interface
    std::visit([this, &U](const auto& eos) {
        std::visit([this, &eos, &U](const auto& flux_scheme) {
            const int first = mesh_.first_interior();
            const int last = mesh_.last_interior();

            // Loop over interfaces (from first interior left face to last interior right face)
            for (int i = first - 1; i <= last; ++i) {
                ConservativeVars U_L, U_R;

                if (order_ >= 2) {
                    // MUSCL reconstruction
                    auto [W_L, W_R] = MUSCLReconstruction::reconstruct(
                        std::span<const PrimitiveVars>(W_), i, limiter_);
                    U_L = eos.to_conservative(W_L);
                    U_R = eos.to_conservative(W_R);
                } else {
                    // First order: piecewise constant
                    U_L = U[static_cast<std::size_t>(i)];
                    U_R = U[static_cast<std::size_t>(i + 1)];
                }

                // Compute numerical flux
                fluxes_[static_cast<std::size_t>(i + 1)] = flux_scheme(U_L, U_R, eos);
            }
        }, flux_);
    }, eos_);

    // Compute dU/dt = -dF/dx = -(F_{i+1/2} - F_{i-1/2}) / dx
    const Real inv_dx = Real{1} / mesh_.dx();
    const int first = mesh_.first_interior();
    const int last = mesh_.last_interior();

    // Zero out dU
    for (std::size_t i = 0; i < dU.size(); ++i) {
        dU[i] = ConservativeVars{};
    }

    // Interior cells only
    for (int i = first; i <= last; ++i) {
        const auto& F_right = fluxes_[static_cast<std::size_t>(i + 1)];
        const auto& F_left = fluxes_[static_cast<std::size_t>(i)];
        dU[static_cast<std::size_t>(i)] = (F_left - F_right) * inv_dx;
    }
}

PrimitiveArray Solver::to_primitive() const {
    PrimitiveArray W(U_.size());
    std::visit([&W, this](const auto& eos) {
        for (std::size_t i = 0; i < U_.size(); ++i) {
            W[i] = eos.to_primitive(U_[i]);
        }
    }, eos_);
    return W;
}

void Solver::run() {
    const Real t_final = config_.time.final_time;
    int step = 0;

    std::println("Starting simulation: {}", config_.simulation.test_name);
    std::println("  Domain: [{}, {}], Cells: {}", mesh_.xmin(), mesh_.xmax(), mesh_.num_cells());
    std::println("  Final time: {}, CFL: {}", t_final, config_.time.cfl);
    std::println("  Order: {}", order_);

    // Start timing
    const auto start_time = std::chrono::high_resolution_clock::now();

    // Define RHS function for time integrator
    auto rhs_func = [this](std::span<const ConservativeVars> U_in, std::span<ConservativeVars> dU_out) {
        // Need a mutable copy for boundary application
        ConservativeArray U_temp(U_in.begin(), U_in.end());
        apply_left_boundary(bc_left_, U_temp, mesh_);
        apply_right_boundary(bc_right_, U_temp, mesh_);
        compute_rhs(U_temp, dU_out);
    };

    while (time_ < t_final) {
        // Compute stable timestep
        Real dt = compute_dt();

        // Adjust final step to hit t_final exactly
        if (time_ + dt > t_final) {
            dt = t_final - time_;
        }

        // Advance solution
        advance(time_integrator_, U_, dt, rhs_func);

        // Apply boundary conditions
        apply_boundaries();

        time_ += dt;
        ++step;

        // Progress output every 100 steps
        if (step % 100 == 0) {
            std::println("  Step {:6d}, t = {:.6f}, dt = {:.6e}", step, time_, dt);
        }
    }

    // End timing and compute performance metrics
    const auto end_time = std::chrono::high_resolution_clock::now();
    const auto elapsed = std::chrono::duration<double>(end_time - start_time);
    const double wall_time = elapsed.count();
    const double cells_per_sec = static_cast<double>(step) * static_cast<double>(mesh_.num_cells()) / wall_time;
    const double steps_per_sec = static_cast<double>(step) / wall_time;

    std::println("Simulation complete: {} steps, final time = {:.6f}", step, time_);
    std::println("Performance:");
    std::println("  Wall time:    {:.4f} s", wall_time);
    std::println("  Steps/sec:    {:.2f}", steps_per_sec);
    std::println("  Mcells/sec:   {:.2f}", cells_per_sec / 1.0e6);
}

}  // namespace euler1d
