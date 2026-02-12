/**
 * @file time_integrator.hpp
 * @brief Time integration schemes for the 1D Euler solver
 */

#ifndef EULER1D_TIME_TIME_INTEGRATOR_HPP
#define EULER1D_TIME_TIME_INTEGRATOR_HPP

#include "../core/types.hpp"
#include <functional>
#include <span>
#include <variant>

namespace euler1d {

/// Type alias for the RHS function: computes dU/dt given U
using RhsFunction = std::function<void(std::span<const ConservativeVars>, std::span<ConservativeVars>)>;

// =============================================================================
// Explicit (Forward) Euler
// =============================================================================

/**
 * @brief Forward Euler time integrator (first order)
 *
 * U^{n+1} = U^n + dt * L(U^n)
 */
struct ExplicitEuler {
    void advance(std::span<ConservativeVars> U, Real dt, const RhsFunction& rhs) const {
        const std::size_t n = U.size();
        ConservativeArray dU(n);

        // Compute RHS
        rhs(U, dU);

        // Update solution
        for (std::size_t i = 0; i < n; ++i) {
            U[i] += dt * dU[i];
        }
    }
};

// =============================================================================
// Strong Stability Preserving RK3 (SSPRK3)
// =============================================================================

/**
 * @brief SSPRK3 time integrator (third order)
 *
 * Three-stage strong stability preserving Runge-Kutta method.
 * Maintains TVD property when spatial discretization is TVD.
 *
 * U^(1) = U^n + dt * L(U^n)
 * U^(2) = 3/4 * U^n + 1/4 * U^(1) + 1/4 * dt * L(U^(1))
 * U^(n+1) = 1/3 * U^n + 2/3 * U^(2) + 2/3 * dt * L(U^(2))
 */
struct SSPRK3 {
    void advance(std::span<ConservativeVars> U, Real dt, const RhsFunction& rhs) const {
        const std::size_t n = U.size();

        // Storage for stages
        ConservativeArray U_n(n);     // U^n
        ConservativeArray U_1(n);     // U^(1)
        ConservativeArray U_2(n);     // U^(2)
        ConservativeArray dU(n);      // RHS

        // Save initial state
        std::copy(U.begin(), U.end(), U_n.begin());

        // Stage 1: U^(1) = U^n + dt * L(U^n)
        rhs(U, dU);
        for (std::size_t i = 0; i < n; ++i) {
            U_1[i] = U_n[i] + dt * dU[i];
        }

        // Stage 2: U^(2) = 3/4 * U^n + 1/4 * U^(1) + 1/4 * dt * L(U^(1))
        rhs(U_1, dU);
        for (std::size_t i = 0; i < n; ++i) {
            U_2[i] = Real{0.75} * U_n[i] + Real{0.25} * U_1[i] + Real{0.25} * dt * dU[i];
        }

        // Stage 3: U^(n+1) = 1/3 * U^n + 2/3 * U^(2) + 2/3 * dt * L(U^(2))
        rhs(U_2, dU);
        for (std::size_t i = 0; i < n; ++i) {
            U[i] = Real{1.0 / 3.0} * U_n[i] + Real{2.0 / 3.0} * U_2[i] + Real{2.0 / 3.0} * dt * dU[i];
        }
    }
};

// =============================================================================
// Time Integrator Variant for runtime selection
// =============================================================================

/// Variant holding all supported time integrators
using TimeIntegratorVariant = std::variant<ExplicitEuler, SSPRK3>;

/// Advance solution by one timestep
inline void advance(const TimeIntegratorVariant& integrator, std::span<ConservativeVars> U,
                   Real dt, const RhsFunction& rhs) {
    std::visit([&](const auto& integ) { integ.advance(U, dt, rhs); }, integrator);
}

}  // namespace euler1d

#endif  // EULER1D_TIME_TIME_INTEGRATOR_HPP
