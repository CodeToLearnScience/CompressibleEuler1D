/**
 * @file types.hpp
 * @brief Core type definitions for the 1D Euler solver
 *
 * Defines the fundamental types used throughout the solver:
 * - Real: compile-time selectable precision (float/double)
 * - ConservativeVars: conserved variables (rho, rho*u, E)
 * - PrimitiveVars: primitive variables (rho, u, p)
 */

#ifndef EULER1D_CORE_TYPES_HPP
#define EULER1D_CORE_TYPES_HPP

#include <array>
#include <cmath>
#include <concepts>
#include <vector>

namespace euler1d {

// =============================================================================
// Compile-time precision selection
// =============================================================================

// Use token pasting to check precision
#define EULER1D_PRECISION_float 1
#define EULER1D_PRECISION_double 2

#define EULER1D_CONCAT_(a, b) a ## b
#define EULER1D_CONCAT(a, b) EULER1D_CONCAT_(a, b)

#if defined(EULER1D_PRECISION) && EULER1D_CONCAT(EULER1D_PRECISION_, EULER1D_PRECISION) == 1
    using Real = float;
#else
    using Real = double;
#endif

// =============================================================================
// Conservative Variables: (rho, rho*u, E)
// =============================================================================

/**
 * @brief Conservative variables for the 1D Euler equations
 *
 * The 1D Euler equations in conservative form:
 *   ∂U/∂t + ∂F(U)/∂x = 0
 *
 * where U = (ρ, ρu, E)^T
 */
struct ConservativeVars {
    Real rho;    ///< Density
    Real rho_u;  ///< Momentum (density * velocity)
    Real E;      ///< Total energy per unit volume

    /// Default constructor (zero initialization)
    constexpr ConservativeVars() noexcept : rho{0}, rho_u{0}, E{0} {}

    /// Value constructor
    constexpr ConservativeVars(Real rho_, Real rho_u_, Real E_) noexcept
        : rho{rho_}, rho_u{rho_u_}, E{E_} {}

    /// Arithmetic operations
    constexpr ConservativeVars operator+(const ConservativeVars& other) const noexcept {
        return {rho + other.rho, rho_u + other.rho_u, E + other.E};
    }

    constexpr ConservativeVars operator-(const ConservativeVars& other) const noexcept {
        return {rho - other.rho, rho_u - other.rho_u, E - other.E};
    }

    constexpr ConservativeVars operator*(Real scalar) const noexcept {
        return {rho * scalar, rho_u * scalar, E * scalar};
    }

    constexpr ConservativeVars operator/(Real scalar) const noexcept {
        return {rho / scalar, rho_u / scalar, E / scalar};
    }

    constexpr ConservativeVars& operator+=(const ConservativeVars& other) noexcept {
        rho += other.rho;
        rho_u += other.rho_u;
        E += other.E;
        return *this;
    }

    constexpr ConservativeVars& operator-=(const ConservativeVars& other) noexcept {
        rho -= other.rho;
        rho_u -= other.rho_u;
        E -= other.E;
        return *this;
    }

    constexpr ConservativeVars& operator*=(Real scalar) noexcept {
        rho *= scalar;
        rho_u *= scalar;
        E *= scalar;
        return *this;
    }

    /// Access by index (0=rho, 1=rho_u, 2=E)
    constexpr Real& operator[](std::size_t i) noexcept {
        switch (i) {
            case 0: return rho;
            case 1: return rho_u;
            default: return E;
        }
    }

    constexpr const Real& operator[](std::size_t i) const noexcept {
        switch (i) {
            case 0: return rho;
            case 1: return rho_u;
            default: return E;
        }
    }

    /// Number of components
    static constexpr std::size_t size() noexcept { return 3; }
};

/// Scalar multiplication (scalar * vars)
constexpr ConservativeVars operator*(Real scalar, const ConservativeVars& vars) noexcept {
    return vars * scalar;
}

// =============================================================================
// Primitive Variables: (rho, u, p)
// =============================================================================

/**
 * @brief Primitive variables for the 1D Euler equations
 *
 * Primitive form: (ρ, u, p)^T
 * More intuitive and often used for reconstruction
 */
struct PrimitiveVars {
    Real rho;  ///< Density
    Real u;    ///< Velocity
    Real p;    ///< Pressure

    /// Default constructor (zero initialization)
    constexpr PrimitiveVars() noexcept : rho{0}, u{0}, p{0} {}

    /// Value constructor
    constexpr PrimitiveVars(Real rho_, Real u_, Real p_) noexcept
        : rho{rho_}, u{u_}, p{p_} {}

    /// Arithmetic operations
    constexpr PrimitiveVars operator+(const PrimitiveVars& other) const noexcept {
        return {rho + other.rho, u + other.u, p + other.p};
    }

    constexpr PrimitiveVars operator-(const PrimitiveVars& other) const noexcept {
        return {rho - other.rho, u - other.u, p - other.p};
    }

    constexpr PrimitiveVars operator*(Real scalar) const noexcept {
        return {rho * scalar, u * scalar, p * scalar};
    }

    constexpr PrimitiveVars operator/(Real scalar) const noexcept {
        return {rho / scalar, u / scalar, p / scalar};
    }

    constexpr PrimitiveVars& operator+=(const PrimitiveVars& other) noexcept {
        rho += other.rho;
        u += other.u;
        p += other.p;
        return *this;
    }

    constexpr PrimitiveVars& operator-=(const PrimitiveVars& other) noexcept {
        rho -= other.rho;
        u -= other.u;
        p -= other.p;
        return *this;
    }

    constexpr PrimitiveVars& operator*=(Real scalar) noexcept {
        rho *= scalar;
        u *= scalar;
        p *= scalar;
        return *this;
    }

    /// Access by index (0=rho, 1=u, 2=p)
    constexpr Real& operator[](std::size_t i) noexcept {
        switch (i) {
            case 0: return rho;
            case 1: return u;
            default: return p;
        }
    }

    constexpr const Real& operator[](std::size_t i) const noexcept {
        switch (i) {
            case 0: return rho;
            case 1: return u;
            default: return p;
        }
    }

    /// Number of components
    static constexpr std::size_t size() noexcept { return 3; }
};

/// Scalar multiplication (scalar * vars)
constexpr PrimitiveVars operator*(Real scalar, const PrimitiveVars& vars) noexcept {
    return vars * scalar;
}

// =============================================================================
// Type aliases for solution arrays
// =============================================================================

/// Array of conservative variables (one per cell including ghosts)
using ConservativeArray = std::vector<ConservativeVars>;

/// Array of primitive variables
using PrimitiveArray = std::vector<PrimitiveVars>;

}  // namespace euler1d

#endif  // EULER1D_CORE_TYPES_HPP
