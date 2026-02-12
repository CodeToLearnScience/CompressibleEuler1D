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

/// Helper to select precision-appropriate constants
template<typename T>
struct PrecisionTraits {
    static constexpr T epsilon = T(1.0e-14);
    static constexpr T min_value = T(1.0e-14);
};

template<>
struct PrecisionTraits<float> {
    static constexpr float epsilon = 1.0e-7f;
    static constexpr float min_value = 1.0e-10f;
};

/// Small number to avoid division by zero
inline constexpr Real epsilon = PrecisionTraits<Real>::epsilon;

/// Default gamma for ideal gas (air at standard conditions)
inline constexpr Real default_gamma = static_cast<Real>(1.4);

/// Minimum allowed density
inline constexpr Real min_density = PrecisionTraits<Real>::min_value;

/// Minimum allowed pressure
inline constexpr Real min_pressure = PrecisionTraits<Real>::min_value;

}  // namespace constants

}  // namespace euler1d

#endif  // EULER1D_CORE_CONSTANTS_HPP
