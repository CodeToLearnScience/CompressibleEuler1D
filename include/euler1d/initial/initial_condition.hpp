/**
 * @file initial_condition.hpp
 * @brief Initial condition generators for the 1D Euler solver
 */

#ifndef EULER1D_INITIAL_INITIAL_CONDITION_HPP
#define EULER1D_INITIAL_INITIAL_CONDITION_HPP

#include "../core/types.hpp"
#include "../core/constants.hpp"
#include "../config/config_types.hpp"
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"
#include <cmath>
#include <span>
#include <variant>

namespace euler1d {

// =============================================================================
// Piecewise Constant Initial Condition
// =============================================================================

/**
 * @brief Piecewise constant initial condition
 *
 * Initializes solution with constant values in specified regions.
 */
struct PiecewiseConstantIC {
    std::vector<Region> regions;

    template <typename Eos>
    void apply(std::span<ConservativeVars> U, const Mesh1D& mesh, const Eos& eos) const {
        for (int i = 0; i < mesh.total_cells(); ++i) {
            const Real x = mesh.x(i);

            // Find which region this cell belongs to
            PrimitiveVars W{Real{1}, Real{0}, Real{1}};  // Default
            for (const auto& region : regions) {
                if (x >= region.x_left && x < region.x_right) {
                    W = PrimitiveVars{region.rho, region.u, region.p};
                    break;
                }
            }

            U[static_cast<std::size_t>(i)] = eos.to_conservative(W);
        }
    }
};

// =============================================================================
// Shock-Entropy Interaction Initial Condition
// =============================================================================

/**
 * @brief Shock-entropy wave interaction initial condition
 *
 * Left of discontinuity: constant state
 * Right of discontinuity: sinusoidal density perturbation
 */
struct ShockEntropyInteractionIC {
    Real discontinuity_position;
    ConstantState left_state;
    SinusoidalState right_state;

    template <typename Eos>
    void apply(std::span<ConservativeVars> U, const Mesh1D& mesh, const Eos& eos) const {
        for (int i = 0; i < mesh.total_cells(); ++i) {
            const Real x = mesh.x(i);
            PrimitiveVars W;

            if (x < discontinuity_position) {
                // Constant left state
                W = PrimitiveVars{left_state.rho, left_state.u, left_state.p};
            } else {
                // Sinusoidal right state
                Real arg = right_state.rho_frequency * x;
                if (right_state.use_pi) {
                    arg *= constants::pi;
                }
                const Real rho = right_state.rho_base + right_state.rho_amplitude * std::sin(arg);
                W = PrimitiveVars{rho, right_state.u, right_state.p};
            }

            U[static_cast<std::size_t>(i)] = eos.to_conservative(W);
        }
    }
};

// =============================================================================
// Initial Condition Variant for runtime selection
// =============================================================================

/// Variant holding all supported initial conditions
using InitialConditionVariant = std::variant<PiecewiseConstantIC, ShockEntropyInteractionIC>;

/// Apply initial condition with any EOS
template <typename Eos>
void apply_initial_condition(const InitialConditionVariant& ic, std::span<ConservativeVars> U,
                              const Mesh1D& mesh, const Eos& eos) {
    std::visit([&](const auto& condition) { condition.apply(U, mesh, eos); }, ic);
}

/// Create initial condition from config
InitialConditionVariant create_initial_condition(const InitialConditionConfig& config);

}  // namespace euler1d

#endif  // EULER1D_INITIAL_INITIAL_CONDITION_HPP
