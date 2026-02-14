/**
 * @file flux.hpp
 * @brief Numerical flux schemes for the 1D Euler equations
 *
 * All flux schemes compute the numerical flux at a cell interface
 * given left and right states.
 */

#ifndef EULER1D_FLUX_FLUX_HPP
#define EULER1D_FLUX_FLUX_HPP

#include "../core/types.hpp"
#include "../eos/eos.hpp"
#include <algorithm>
#include <cmath>
#include <variant>

namespace euler1d {

// =============================================================================
// Local Lax-Friedrichs (LLF) Flux
// =============================================================================

/**
 * @brief Local Lax-Friedrichs (Rusanov-type) flux
 *
 * F_{i+1/2} = 0.5 * (F_L + F_R) - 0.5 * λ_max * (U_R - U_L)
 * where λ_max = max(|u_L| + c_L, |u_R| + c_R)
 */
struct LLFFlux {
    template <typename Eos>
    [[nodiscard]] ConservativeVars operator()(
        const ConservativeVars& U_L,
        const ConservativeVars& U_R,
        const Eos& eos) const noexcept {

        // Compute physical fluxes
        const auto F_L = eos.flux(U_L);
        const auto F_R = eos.flux(U_R);

        // Compute wave speeds
        const Real u_L = U_L.rho_u / U_L.rho;
        const Real u_R = U_R.rho_u / U_R.rho;
        const Real c_L = eos.sound_speed(U_L);
        const Real c_R = eos.sound_speed(U_R);

        // Maximum wave speed
        const Real lambda_max = std::max(std::abs(u_L) + c_L, std::abs(u_R) + c_R);

        // LLF flux
        return Real{0.5} * (F_L + F_R) - Real{0.5} * lambda_max * (U_R - U_L);
    }
};

// =============================================================================
// Rusanov Flux (identical to LLF)
// =============================================================================

/**
 * @brief Rusanov flux (alias for LLF)
 */
struct RusanovFlux {
    template <typename Eos>
    [[nodiscard]] ConservativeVars operator()(
        const ConservativeVars& U_L,
        const ConservativeVars& U_R,
        const Eos& eos) const noexcept {
        return LLFFlux{}(U_L, U_R, eos);
    }
};

// =============================================================================
// HLL Flux (Harten-Lax-van Leer)
// =============================================================================

/**
 * @brief HLL flux with Davis wave speed estimates
 *
 * Two-wave approximate Riemann solver using fastest left and right
 * wave speeds.
 */
struct HLLFlux {
    template <typename Eos>
    [[nodiscard]] ConservativeVars operator()(
        const ConservativeVars& U_L,
        const ConservativeVars& U_R,
        const Eos& eos) const noexcept {

        // Left state
        const Real rho_L = U_L.rho;
        const Real u_L = U_L.rho_u / rho_L;
        const Real p_L = eos.pressure(U_L);
        const Real c_L = eos.sound_speed(rho_L, p_L);

        // Right state
        const Real rho_R = U_R.rho;
        const Real u_R = U_R.rho_u / rho_R;
        const Real p_R = eos.pressure(U_R);
        const Real c_R = eos.sound_speed(rho_R, p_R);

        // Davis wave speed estimates
        const Real S_L = std::min(u_L - c_L, u_R - c_R);
        const Real S_R = std::max(u_L + c_L, u_R + c_R);

        // Physical fluxes
        const auto F_L = eos.flux(U_L);
        const auto F_R = eos.flux(U_R);

        // HLL flux
        if (S_L >= Real{0}) {
            return F_L;
        } else if (S_R <= Real{0}) {
            return F_R;
        } else {
            return (S_R * F_L - S_L * F_R + S_L * S_R * (U_R - U_L)) / (S_R - S_L);
        }
    }
};

// =============================================================================
// HLLC Flux (HLL with Contact restoration)
// =============================================================================

/**
 * @brief HLLC flux - three-wave approximate Riemann solver
 *
 * Restores the contact discontinuity missing in HLL, providing
 * better resolution of contact waves and shear layers.
 */
struct HLLCFlux {
    template <typename Eos>
    [[nodiscard]] ConservativeVars operator()(
        const ConservativeVars& U_L,
        const ConservativeVars& U_R,
        const Eos& eos) const noexcept {

        // Left state
        const Real rho_L = U_L.rho;
        const Real u_L = U_L.rho_u / rho_L;
        const Real p_L = eos.pressure(U_L);
        const Real c_L = eos.sound_speed(rho_L, p_L);
        const Real E_L = U_L.E;

        // Right state
        const Real rho_R = U_R.rho;
        const Real u_R = U_R.rho_u / rho_R;
        const Real p_R = eos.pressure(U_R);
        const Real c_R = eos.sound_speed(rho_R, p_R);
        const Real E_R = U_R.E;

        // Wave speed estimates (Davis estimates)
        const Real S_L = std::min(u_L - c_L, u_R - c_R);
        const Real S_R = std::max(u_L + c_L, u_R + c_R);

        // Contact wave speed
        const Real S_star = (p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)) /
                            (rho_L * (S_L - u_L) - rho_R * (S_R - u_R));

        // Physical fluxes
        const auto F_L = eos.flux(U_L);
        const auto F_R = eos.flux(U_R);

        if (S_L >= Real{0}) {
            return F_L;
        } else if (S_R <= Real{0}) {
            return F_R;
        } else if (S_star >= Real{0}) {
            // Left star state
            const Real coeff = rho_L * (S_L - u_L) / (S_L - S_star);
            const ConservativeVars U_star_L{
                coeff,
                coeff * S_star,
                coeff * (E_L / rho_L + (S_star - u_L) * (S_star + p_L / (rho_L * (S_L - u_L))))
            };
            return F_L + S_L * (U_star_L - U_L);
        } else {
            // Right star state
            const Real coeff = rho_R * (S_R - u_R) / (S_R - S_star);
            const ConservativeVars U_star_R{
                coeff,
                coeff * S_star,
                coeff * (E_R / rho_R + (S_star - u_R) * (S_star + p_R / (rho_R * (S_R - u_R))))
            };
            return F_R + S_R * (U_star_R - U_R);
        }
    }
};

// =============================================================================
// MOVERS Flux (Exact shock and contact wave speed estimates)
// =============================================================================

/**
 * @brief MOVERS flux - exact shock and contact wave speed estimates
 *
 * Steady contacts and shocks are captured exactly 
 */
struct MoversLEFlux {
    template <typename Eos>
    [[nodiscard]] ConservativeVars operator()(
        const ConservativeVars& U_L,
        const ConservativeVars& U_R,
        const Eos& eos) const noexcept {

        // Compute physical fluxes
        const auto F_L = eos.flux(U_L);
        const auto F_R = eos.flux(U_R);

        // Compute wave speeds
        const Real u_L = U_L.rho_u / U_L.rho;
        const Real u_R = U_R.rho_u / U_R.rho;
        const Real c_L = eos.sound_speed(U_L);
        const Real c_R = eos.sound_speed(U_R);

        // Maximum Minimum wave speed
        const auto lambda_max_min = max_min_eig_value(u_L, u_R, c_L, c_R);

        // MoversLE flux (component-wise dissipation)
        auto flux_component = [this, &lambda_max_min](Real flux_R, Real flux_L, Real U_R_var, Real U_L_var) {
            const Real diss = compute_dissipation(flux_R, flux_L, U_R_var, U_L_var, lambda_max_min.first, lambda_max_min.second);
            return Real{0.5} * (flux_L + flux_R) - Real{0.5} * diss * (U_R_var - U_L_var);
        };

        return {
            flux_component(F_R.rho, F_L.rho, U_R.rho, U_L.rho),
            flux_component(F_R.rho_u, F_L.rho_u, U_R.rho_u, U_L.rho_u),
            flux_component(F_R.E, F_L.E, U_R.E, U_L.E)
        };
    }

    auto max_min_eig_value(Real ul, Real ur, Real cl, Real cr) const noexcept {
        const Real L1r = std::abs(ur + cr), L1l = std::abs(ul + cl);
        const Real L2r = std::abs(ur),     L2l = std::abs(ul);
        const Real L3r = std::abs(ur - cr), L3l = std::abs(ul - cl);
        
        auto max_eig_va = std::max({std::max({L1r, L2r, L3r}), std::max({L1l, L2l, L3l})});    
        auto min_eig_val =  std::min({std::min({L1r, L2r, L3r}), std::min({L1l, L2l, L3l})});
        return std::pair<Real, Real>{max_eig_va, min_eig_val};
    }

    auto compute_dissipation(Real Fr, Real Fl, Real Ur, Real Ul, Real L_max, Real L_min) const noexcept {
        const Real epsilon = Real{1e-6};
        Real S{0.0};
        if (std::abs(Fr - Fl) < epsilon) return Real{0};
        if (std::abs(Ur - Ul) < epsilon) return L_min;
        
        if (std::abs(Ur - Ul) > epsilon && std::abs(Fr - Fl) > epsilon) {
            S = std::abs((Fr - Fl) / (Ur - Ul));
        } else {
            S = L_min;
        }
        
        if (S < epsilon) return Real{0};
        if (S >= L_max) return L_max;
        if (S <= L_min) return L_min;
        return S;
    }
};

// =============================================================================
// Flux Variant for runtime selection
// =============================================================================

/// Variant holding all supported numerical flux schemes
using FluxVariant = std::variant<LLFFlux, RusanovFlux, HLLFlux, HLLCFlux, MoversLEFlux>;

/// Compute numerical flux using any flux scheme
template <typename Eos>
[[nodiscard]] inline ConservativeVars compute_flux(
    const FluxVariant& flux,
    const ConservativeVars& U_L,
    const ConservativeVars& U_R,
    const Eos& eos) {
    return std::visit([&](const auto& f) { return f(U_L, U_R, eos); }, flux);
}

}  // namespace euler1d

#endif  // EULER1D_FLUX_FLUX_HPP
