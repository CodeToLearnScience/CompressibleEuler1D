/**
 * @file limiter.hpp
 * @brief Slope limiters for MUSCL reconstruction
 *
 * TVD slope limiters to prevent spurious oscillations near discontinuities.
 */

#ifndef EULER1D_RECONSTRUCTION_LIMITER_HPP
#define EULER1D_RECONSTRUCTION_LIMITER_HPP

#include "../core/types.hpp"
#include <algorithm>
#include <cmath>
#include <variant>

namespace euler1d {

// =============================================================================
// No Limiter (for first-order or testing)
// =============================================================================

/**
 * @brief No limiting - returns 0 (first order) or 1 (central diff)
 */
struct NoLimiter {
    [[nodiscard]] constexpr Real operator()(Real /*r*/) const noexcept {
        return Real{0};  // First order
    }
};

// =============================================================================
// Minmod Limiter
// =============================================================================

/**
 * @brief Minmod limiter - most diffusive TVD limiter
 *
 * φ(r) = max(0, min(1, r))
 */
struct MinmodLimiter {
    [[nodiscard]] constexpr Real operator()(Real r) const noexcept {
        return std::max(Real{0}, std::min(Real{1}, r));
    }
};

// =============================================================================
// Van Leer Limiter
// =============================================================================

/**
 * @brief Van Leer limiter - smooth TVD limiter
 *
 * φ(r) = (r + |r|) / (1 + |r|)
 */
struct VanLeerLimiter {
    [[nodiscard]] Real operator()(Real r) const noexcept {
        return (r + std::abs(r)) / (Real{1} + std::abs(r));
    }
};

// =============================================================================
// Superbee Limiter
// =============================================================================

/**
 * @brief Superbee limiter - least diffusive TVD limiter
 *
 * φ(r) = max(0, min(2r, 1), min(r, 2))
 */
struct SuperbeeLimiter {
    [[nodiscard]] constexpr Real operator()(Real r) const noexcept {
        return std::max({Real{0}, std::min(Real{2} * r, Real{1}), std::min(r, Real{2})});
    }
};

// =============================================================================
// MC (Monotonized Central) Limiter
// =============================================================================

/**
 * @brief MC (Monotonized Central) limiter
 *
 * φ(r) = max(0, min(2r, (1+r)/2, 2))
 */
struct MCLimiter {
    [[nodiscard]] constexpr Real operator()(Real r) const noexcept {
        return std::max(Real{0}, std::min({Real{2} * r, (Real{1} + r) / Real{2}, Real{2}}));
    }
};

// =============================================================================
// Limiter Variant for runtime selection
// =============================================================================

/// Variant holding all supported limiters
using LimiterVariant = std::variant<NoLimiter, MinmodLimiter, VanLeerLimiter, SuperbeeLimiter, MCLimiter>;

/// Apply limiter function
[[nodiscard]] inline Real apply_limiter(const LimiterVariant& limiter, Real r) {
    return std::visit([r](const auto& lim) { return lim(r); }, limiter);
}

}  // namespace euler1d

#endif  // EULER1D_RECONSTRUCTION_LIMITER_HPP
