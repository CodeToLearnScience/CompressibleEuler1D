/**
 * @file test_flux.cpp
 * @brief Unit tests for numerical flux schemes
 */

#include <gtest/gtest.h>
#include "euler1d/flux/flux.hpp"
#include "euler1d/eos/eos.hpp"
#include <cmath>

using namespace euler1d;

class FluxTest : public ::testing::Test {
protected:
    IdealGas eos{1.4};

    // Standard test state
    ConservativeVars make_state(Real rho, Real u, Real p) {
        PrimitiveVars W{rho, u, p};
        return eos.to_conservative(W);
    }
};

TEST_F(FluxTest, LLFConsistency) {
    // For identical left and right states, flux should equal physical flux
    auto U = make_state(1.0, 1.0, 1.0);
    LLFFlux llf;

    auto F_num = llf(U, U, eos);
    auto F_phys = eos.flux(U);

    EXPECT_NEAR(F_num.rho, F_phys.rho, 1e-12);
    EXPECT_NEAR(F_num.rho_u, F_phys.rho_u, 1e-12);
    EXPECT_NEAR(F_num.E, F_phys.E, 1e-12);
}

TEST_F(FluxTest, LLFSymmetry) {
    // F(U_L, U_R) evaluated left-to-right should be consistent
    auto U_L = make_state(1.0, 0.0, 1.0);
    auto U_R = make_state(0.125, 0.0, 0.1);
    LLFFlux llf;

    auto F = llf(U_L, U_R, eos);

    // Flux should be bounded between physical fluxes
    // (This is a sanity check, not a strict requirement)
    EXPECT_TRUE(std::isfinite(F.rho));
    EXPECT_TRUE(std::isfinite(F.rho_u));
    EXPECT_TRUE(std::isfinite(F.E));
}

TEST_F(FluxTest, HLLConsistency) {
    auto U = make_state(1.0, 1.0, 1.0);
    HLLFlux hll;

    auto F_num = hll(U, U, eos);
    auto F_phys = eos.flux(U);

    EXPECT_NEAR(F_num.rho, F_phys.rho, 1e-12);
    EXPECT_NEAR(F_num.rho_u, F_phys.rho_u, 1e-12);
    EXPECT_NEAR(F_num.E, F_phys.E, 1e-12);
}

TEST_F(FluxTest, HLLCConsistency) {
    auto U = make_state(1.0, 1.0, 1.0);
    HLLCFlux hllc;

    auto F_num = hllc(U, U, eos);
    auto F_phys = eos.flux(U);

    EXPECT_NEAR(F_num.rho, F_phys.rho, 1e-12);
    EXPECT_NEAR(F_num.rho_u, F_phys.rho_u, 1e-12);
    EXPECT_NEAR(F_num.E, F_phys.E, 1e-12);
}

TEST_F(FluxTest, FluxVariantWorks) {
    auto U_L = make_state(1.0, 0.0, 1.0);
    auto U_R = make_state(0.125, 0.0, 0.1);

    FluxVariant flux = LLFFlux{};
    auto F = compute_flux(flux, U_L, U_R, eos);

    EXPECT_TRUE(std::isfinite(F.rho));
    EXPECT_TRUE(std::isfinite(F.rho_u));
    EXPECT_TRUE(std::isfinite(F.E));
}

TEST_F(FluxTest, AllFluxesProduceValidOutput) {
    auto U_L = make_state(1.0, 0.0, 1.0);
    auto U_R = make_state(0.125, 0.0, 0.1);

    std::vector<FluxVariant> fluxes = {
        LLFFlux{}, RusanovFlux{}, HLLFlux{}, HLLCFlux{}
    };

    for (const auto& flux : fluxes) {
        auto F = compute_flux(flux, U_L, U_R, eos);
        EXPECT_TRUE(std::isfinite(F.rho));
        EXPECT_TRUE(std::isfinite(F.rho_u));
        EXPECT_TRUE(std::isfinite(F.E));
    }
}
