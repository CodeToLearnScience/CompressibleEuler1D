/**
 * @file test_reconstruction.cpp
 * @brief Unit tests for MUSCL reconstruction
 */

#include <gtest/gtest.h>
#include "euler1d/reconstruction/muscl.hpp"

using namespace euler1d;

TEST(ReconstructionTest, FirstOrderReturnsOriginalStates) {
    PrimitiveArray W = {
        {1.0, 0.0, 1.0},
        {1.0, 0.0, 1.0},
        {0.5, 0.0, 0.5},
        {0.5, 0.0, 0.5}
    };

    auto [W_L, W_R] = FirstOrderReconstruction::reconstruct(W, 1);

    EXPECT_DOUBLE_EQ(W_L.rho, W[1].rho);
    EXPECT_DOUBLE_EQ(W_R.rho, W[2].rho);
}

TEST(ReconstructionTest, MUSCLProducesValidStates) {
    PrimitiveArray W = {
        {1.0, 1.0, 1.0},
        {1.0, 1.0, 1.0},
        {1.0, 1.0, 1.0},
        {0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5}
    };

    VanLeerLimiter limiter;
    auto [W_L, W_R] = MUSCLReconstruction::reconstruct(
        std::span<const PrimitiveVars>(W), 2, limiter);

    // Reconstructed values should be bounded
    EXPECT_GT(W_L.rho, 0.0);
    EXPECT_GT(W_R.rho, 0.0);
    EXPECT_TRUE(std::isfinite(W_L.rho));
    EXPECT_TRUE(std::isfinite(W_R.rho));
}

TEST(ReconstructionTest, UniformFieldUnchanged) {
    // For uniform field, reconstruction should return same values
    PrimitiveArray W(6, {1.0, 2.0, 3.0});

    MinmodLimiter limiter;
    auto [W_L, W_R] = MUSCLReconstruction::reconstruct(
        std::span<const PrimitiveVars>(W), 2, limiter);

    EXPECT_NEAR(W_L.rho, 1.0, 1e-12);
    EXPECT_NEAR(W_L.u, 2.0, 1e-12);
    EXPECT_NEAR(W_L.p, 3.0, 1e-12);
}

TEST(ReconstructionTest, LimiterVariantWorks) {
    PrimitiveArray W = {
        {1.0, 0.0, 1.0},
        {1.0, 0.0, 1.0},
        {1.0, 0.0, 1.0},
        {0.5, 0.0, 0.5},
        {0.5, 0.0, 0.5}
    };

    LimiterVariant limiter = VanLeerLimiter{};
    auto [W_L, W_R] = MUSCLReconstruction::reconstruct(
        std::span<const PrimitiveVars>(W), 2, limiter);

    EXPECT_TRUE(std::isfinite(W_L.rho));
    EXPECT_TRUE(std::isfinite(W_R.rho));
}
