/**
 * @file constants.hpp
 * @brief Physical and numerical constants for the 1D Euler solver
 */

#ifndef EULER1D_CORE_CONSTANTS_HPP
#define EULER1D_CORE_CONSTANTS_HPP

#include "types.hpp"
#include <numbers>

namespace euler1d {

/// Mathematical constants with appropriate precision
namespace constants {

/// Pi
inline constexpr Real pi = std::numbers::pi_v<Real>;

/// Small number to avoid division by zero
inline constexpr Real epsilon = std::is_same_v<Real, float> ? 1.0e-7f : 1.0e-14;

/// Default gamma for ideal gas (air at standard conditions)
inline constexpr Real default_gamma = static_cast<Real>(1.4);

/// Minimum allowed density
inline constexpr Real min_density = std::is_same_v<Real, float> ? 1.0e-10f : 1.0e-14;

/// Minimum allowed pressure
inline constexpr Real min_pressure = std::is_same_v<Real, float> ? 1.0e-10f : 1.0e-14;

}  // namespace constants

}  // namespace euler1d

#endif  // EULER1D_CORE_CONSTANTS_HPP
