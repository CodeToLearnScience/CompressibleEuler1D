/**
 * @file test_limiter.cpp
 * @brief Unit tests for slope limiters
 */

#include <gtest/gtest.h>
#include "euler1d/reconstruction/limiter.hpp"

using namespace euler1d;

TEST(LimiterTest, MinmodBounds) {
    MinmodLimiter lim;

    // phi(r) should be in [0, 1] for all r
    EXPECT_DOUBLE_EQ(lim(-1.0), 0.0);
    EXPECT_DOUBLE_EQ(lim(0.0), 0.0);
    EXPECT_DOUBLE_EQ(lim(0.5), 0.5);
    EXPECT_DOUBLE_EQ(lim(1.0), 1.0);
    EXPECT_DOUBLE_EQ(lim(2.0), 1.0);
}

TEST(LimiterTest, VanLeerSymmetric) {
    VanLeerLimiter lim;

    // Van Leer: symmetric, phi(r) = (r + |r|) / (1 + |r|)
    EXPECT_NEAR(lim(1.0), 1.0, 1e-12);  // phi(1) = 1
    EXPECT_DOUBLE_EQ(lim(0.0), 0.0);
    EXPECT_GT(lim(2.0), 0.0);
    EXPECT_LT(lim(2.0), 2.0);
}

TEST(LimiterTest, SuperbeeLeastDiffusive) {
    SuperbeeLimiter superbee;
    MinmodLimiter minmod;
    VanLeerLimiter vanleer;

    // Superbee should be >= minmod for r > 0
    for (Real r = 0.1; r <= 3.0; r += 0.1) {
        EXPECT_GE(superbee(r), minmod(r) - 1e-10);
    }
}

TEST(LimiterTest, MCLimiter) {
    MCLimiter mc;

    EXPECT_DOUBLE_EQ(mc(0.0), 0.0);
    EXPECT_DOUBLE_EQ(mc(1.0), 1.0);
    EXPECT_LE(mc(2.0), 2.0);  // bounded by 2
}

TEST(LimiterTest, NoLimiterReturnsZero) {
    NoLimiter lim;
    EXPECT_DOUBLE_EQ(lim(0.5), 0.0);
    EXPECT_DOUBLE_EQ(lim(1.0), 0.0);
    EXPECT_DOUBLE_EQ(lim(2.0), 0.0);
}

TEST(LimiterTest, TVDRegion) {
    // TVD region: limiter should satisfy 0 <= phi(r) <= 2r and 0 <= phi(r) <= 2
    std::vector<LimiterVariant> limiters = {
        MinmodLimiter{}, VanLeerLimiter{}, SuperbeeLimiter{}, MCLimiter{}
    };

    for (const auto& lim : limiters) {
        for (Real r = 0.1; r <= 3.0; r += 0.1) {
            Real phi = apply_limiter(lim, r);
            EXPECT_GE(phi, 0.0);
            EXPECT_LE(phi, std::min(2.0 * r, 2.0) + 1e-10);
        }
    }
}
