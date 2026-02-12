/**
 * @file eos.hpp
 * @brief Equation of State definitions for the 1D Euler solver
 *
 * Provides various equations of state, starting with ideal gas.
 * Uses std::variant for runtime polymorphism without virtual dispatch.
 */

#ifndef EULER1D_EOS_EOS_HPP
#define EULER1D_EOS_EOS_HPP

#include "../core/types.hpp"
#include "../core/constants.hpp"
#include <cmath>
#include <variant>

namespace euler1d {

// =============================================================================
// Ideal Gas Equation of State
// =============================================================================

/**
 * @brief Ideal gas equation of state
 *
 * p = (γ - 1) * ρ * e
 * where e is the specific internal energy
 */
struct IdealGas {
    Real gamma;  ///< Ratio of specific heats (Cp/Cv)

    /// Construct with given gamma
    explicit constexpr IdealGas(Real gamma_ = constants::default_gamma) noexcept
        : gamma{gamma_} {}

    /// Compute pressure from conservative variables
    [[nodiscard]] constexpr Real pressure(const ConservativeVars& U) const noexcept {
        const Real rho = U.rho;
        const Real u = U.rho_u / rho;
        const Real kinetic = Real{0.5} * rho * u * u;
        const Real internal = U.E - kinetic;
        return (gamma - Real{1}) * internal;
    }

    /// Compute pressure from density and internal energy
    [[nodiscard]] constexpr Real pressure(Real rho, Real e_internal) const noexcept {
        return (gamma - Real{1}) * rho * e_internal;
    }

    /// Compute sound speed from density and pressure
    [[nodiscard]] Real sound_speed(Real rho, Real p) const noexcept {
        return std::sqrt(gamma * p / rho);
    }

    /// Compute sound speed from conservative variables
    [[nodiscard]] Real sound_speed(const ConservativeVars& U) const noexcept {
        return sound_speed(U.rho, pressure(U));
    }

    /// Compute specific internal energy from pressure and density
    [[nodiscard]] constexpr Real internal_energy(Real rho, Real p) const noexcept {
        return p / ((gamma - Real{1}) * rho);
    }

    /// Compute total energy from primitive variables
    [[nodiscard]] constexpr Real total_energy(const PrimitiveVars& W) const noexcept {
        const Real e_internal = internal_energy(W.rho, W.p);
        const Real e_kinetic = Real{0.5} * W.u * W.u;
        return W.rho * (e_internal + e_kinetic);
    }

    /// Compute specific enthalpy h = e + p/rho = (E + p)/rho
    [[nodiscard]] constexpr Real enthalpy(const ConservativeVars& U) const noexcept {
        const Real p = pressure(U);
        return (U.E + p) / U.rho;
    }

    /// Compute specific enthalpy from primitive variables
    [[nodiscard]] constexpr Real enthalpy(const PrimitiveVars& W) const noexcept {
        const Real e_int = internal_energy(W.rho, W.p);
        return e_int + Real{0.5} * W.u * W.u + W.p / W.rho;
    }

    /// Convert primitive to conservative variables
    [[nodiscard]] constexpr ConservativeVars to_conservative(const PrimitiveVars& W) const noexcept {
        return ConservativeVars{
            W.rho,
            W.rho * W.u,
            total_energy(W)
        };
    }

    /// Convert conservative to primitive variables
    [[nodiscard]] constexpr PrimitiveVars to_primitive(const ConservativeVars& U) const noexcept {
        const Real rho = U.rho;
        const Real u = U.rho_u / rho;
        const Real p = pressure(U);
        return PrimitiveVars{rho, u, p};
    }

    /// Compute the physical flux F(U)
    [[nodiscard]] constexpr ConservativeVars flux(const ConservativeVars& U) const noexcept {
        const Real rho = U.rho;
        const Real u = U.rho_u / rho;
        const Real p = pressure(U);
        return ConservativeVars{
            U.rho_u,                    // ρu
            U.rho_u * u + p,            // ρu² + p
            (U.E + p) * u               // (E + p)u
        };
    }

    /// Compute the physical flux from primitive variables
    [[nodiscard]] constexpr ConservativeVars flux(const PrimitiveVars& W) const noexcept {
        const Real E = total_energy(W);
        return ConservativeVars{
            W.rho * W.u,                        // ρu
            W.rho * W.u * W.u + W.p,            // ρu² + p
            (E + W.p) * W.u                     // (E + p)u
        };
    }
};

// =============================================================================
// EOS Variant type for runtime selection
// =============================================================================

/// Variant holding all supported equations of state
using EosVariant = std::variant<IdealGas>;

// =============================================================================
// Free functions dispatched via std::visit
// =============================================================================

/// Compute pressure from conservative variables (any EOS)
[[nodiscard]] inline Real pressure(const EosVariant& eos, const ConservativeVars& U) {
    return std::visit([&U](const auto& e) { return e.pressure(U); }, eos);
}

/// Compute sound speed (any EOS)
[[nodiscard]] inline Real sound_speed(const EosVariant& eos, const ConservativeVars& U) {
    return std::visit([&U](const auto& e) { return e.sound_speed(U); }, eos);
}

/// Convert primitive to conservative (any EOS)
[[nodiscard]] inline ConservativeVars to_conservative(const EosVariant& eos, const PrimitiveVars& W) {
    return std::visit([&W](const auto& e) { return e.to_conservative(W); }, eos);
}

/// Convert conservative to primitive (any EOS)
[[nodiscard]] inline PrimitiveVars to_primitive(const EosVariant& eos, const ConservativeVars& U) {
    return std::visit([&U](const auto& e) { return e.to_primitive(U); }, eos);
}

/// Compute physical flux (any EOS)
[[nodiscard]] inline ConservativeVars flux(const EosVariant& eos, const ConservativeVars& U) {
    return std::visit([&U](const auto& e) { return e.flux(U); }, eos);
}

}  // namespace euler1d

#endif  // EULER1D_EOS_EOS_HPP
