/**
 * @file test_eos.cpp
 * @brief Unit tests for equation of state
 */

#include <gtest/gtest.h>
#include "euler1d/eos/eos.hpp"
#include "euler1d/core/constants.hpp"
#include <cmath>

using namespace euler1d;

class IdealGasTest : public ::testing::Test {
protected:
    IdealGas eos{1.4};
};

TEST_F(IdealGasTest, PressureFromPrimitive) {
    // Given primitive state
    PrimitiveVars W{1.0, 0.0, 1.0};  // rho=1, u=0, p=1

    // Convert to conservative
    auto U = eos.to_conservative(W);

    // Recover pressure
    Real p = eos.pressure(U);

    EXPECT_NEAR(p, 1.0, 1e-12);
}

TEST_F(IdealGasTest, SoundSpeed) {
    // For ideal gas: c = sqrt(gamma * p / rho)
    Real rho = 1.0;
    Real p = 1.0;
    Real c = eos.sound_speed(rho, p);

    Real expected = std::sqrt(1.4 * 1.0 / 1.0);  // sqrt(1.4)
    EXPECT_NEAR(c, expected, 1e-12);
}

TEST_F(IdealGasTest, ConservativePrimitiveRoundTrip) {
    // Start with primitive
    PrimitiveVars W_orig{1.225, 100.0, 101325.0};  // air at sea level

    // Convert to conservative
    auto U = eos.to_conservative(W_orig);

    // Convert back to primitive
    auto W = eos.to_primitive(U);

    EXPECT_NEAR(W.rho, W_orig.rho, 1e-10);
    EXPECT_NEAR(W.u, W_orig.u, 1e-10);
    EXPECT_NEAR(W.p, W_orig.p, 1e-6);  // pressure has larger absolute values
}

TEST_F(IdealGasTest, TotalEnergy) {
    // E = rho * (e_internal + 0.5*u^2)
    // e_internal = p / ((gamma-1)*rho)
    PrimitiveVars W{1.0, 10.0, 1.0};

    Real e_int = 1.0 / (0.4 * 1.0);  // p / ((gamma-1)*rho)
    Real e_kin = 0.5 * 10.0 * 10.0;
    Real E_expected = 1.0 * (e_int + e_kin);

    auto U = eos.to_conservative(W);
    EXPECT_NEAR(U.E, E_expected, 1e-10);
}

TEST_F(IdealGasTest, FluxConsistency) {
    // Flux should be consistent with Euler equations
    PrimitiveVars W{1.0, 1.0, 1.0};
    auto U = eos.to_conservative(W);
    auto F = eos.flux(U);

    // F = [rho*u, rho*u^2 + p, (E+p)*u]
    EXPECT_NEAR(F.rho, 1.0, 1e-12);      // rho * u
    EXPECT_NEAR(F.rho_u, 2.0, 1e-12);    // rho*u^2 + p = 1 + 1
    // (E + p) * u
    Real E_plus_p = U.E + 1.0;
    EXPECT_NEAR(F.E, E_plus_p, 1e-10);
}

TEST_F(IdealGasTest, EosVariantWorks) {
    EosVariant eos_var = IdealGas{1.4};

    PrimitiveVars W{1.0, 0.0, 1.0};
    auto U = to_conservative(eos_var, W);

    Real p = pressure(eos_var, U);
    EXPECT_NEAR(p, 1.0, 1e-12);
}
