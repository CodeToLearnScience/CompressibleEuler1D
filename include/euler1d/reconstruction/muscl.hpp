/**
 * @file muscl.hpp
 * @brief MUSCL reconstruction for second-order accuracy
 *
 * Monotone Upstream-centered Schemes for Conservation Laws.
 * Provides piecewise linear reconstruction with slope limiting.
 */

#ifndef EULER1D_RECONSTRUCTION_MUSCL_HPP
#define EULER1D_RECONSTRUCTION_MUSCL_HPP

#include "../core/types.hpp"
#include "../core/constants.hpp"
#include "limiter.hpp"
#include <span>
#include <utility>

namespace euler1d {

/**
 * @brief MUSCL reconstruction on primitive variables
 *
 * Reconstructs left and right states at cell interface i+1/2
 * using piecewise linear reconstruction with slope limiting.
 */
struct MUSCLReconstruction {

    /**
     * @brief Reconstruct primitive variables at interface i+1/2
     *
     * @param W Array of primitive variables (including ghosts)
     * @param i Cell index (left cell of interface)
     * @param limiter Slope limiter to use
     * @return Pair of (W_L, W_R) at interface i+1/2
     */
    template <typename LimiterT>
    [[nodiscard]] static std::pair<PrimitiveVars, PrimitiveVars> reconstruct(
        std::span<const PrimitiveVars> W,
        int i,
        const LimiterT& limiter) {

        // Get stencil values
        const auto& W_im1 = W[static_cast<std::size_t>(i - 1)];  // W_{i-1}
        const auto& W_i   = W[static_cast<std::size_t>(i)];      // W_i
        const auto& W_ip1 = W[static_cast<std::size_t>(i + 1)];  // W_{i+1}
        const auto& W_ip2 = W[static_cast<std::size_t>(i + 2)];  // W_{i+2}

        PrimitiveVars W_L, W_R;

        // Component-wise reconstruction
        for (std::size_t k = 0; k < PrimitiveVars::size(); ++k) {
            // Left state: extrapolate from cell i to right face
            const Real delta_L = W_i[k] - W_im1[k];
            const Real delta_R_left = W_ip1[k] - W_i[k];
            const Real r_L = (std::abs(delta_R_left) > constants::epsilon)
                ? delta_L / delta_R_left
                : Real{0};
            const Real phi_L = limiter(r_L);
            W_L[k] = W_i[k] + Real{0.5} * phi_L * delta_R_left;

            // Right state: extrapolate from cell i+1 to left face
            const Real delta_L_right = W_ip1[k] - W_i[k];
            const Real delta_R = W_ip2[k] - W_ip1[k];
            const Real r_R = (std::abs(delta_L_right) > constants::epsilon)
                ? delta_R / delta_L_right
                : Real{0};
            const Real phi_R = limiter(r_R);
            W_R[k] = W_ip1[k] - Real{0.5} * phi_R * delta_L_right;
        }

        return {W_L, W_R};
    }

    /**
     * @brief Reconstruct with runtime limiter variant
     */
    [[nodiscard]] static std::pair<PrimitiveVars, PrimitiveVars> reconstruct(
        std::span<const PrimitiveVars> W,
        int i,
        const LimiterVariant& limiter) {
        return std::visit([&](const auto& lim) { return reconstruct(W, i, lim); }, limiter);
    }
};

/**
 * @brief First-order reconstruction (no gradients)
 */
struct FirstOrderReconstruction {
    [[nodiscard]] static std::pair<PrimitiveVars, PrimitiveVars> reconstruct(
        std::span<const PrimitiveVars> W,
        int i) {
        return {W[static_cast<std::size_t>(i)], W[static_cast<std::size_t>(i + 1)]};
    }
};

}  // namespace euler1d

#endif  // EULER1D_RECONSTRUCTION_MUSCL_HPP
