/**
 * @file test_time_integrator.cpp
 * @brief Unit tests for time integration schemes
 */

#include <gtest/gtest.h>
#include "euler1d/time/time_integrator.hpp"
#include <cmath>

using namespace euler1d;

class TimeIntegratorTest : public ::testing::Test {
protected:
    // Simple decay ODE: du/dt = -u, solution: u(t) = u0 * exp(-t)
    static void decay_rhs(std::span<const ConservativeVars> U, std::span<ConservativeVars> dU) {
        for (std::size_t i = 0; i < U.size(); ++i) {
            dU[i] = ConservativeVars{-U[i].rho, -U[i].rho_u, -U[i].E};
        }
    }
};

TEST_F(TimeIntegratorTest, ExplicitEulerConverges) {
    ConservativeArray U(1, ConservativeVars{1.0, 0.0, 0.0});
    ExplicitEuler euler;

    Real dt = 0.01;
    int n_steps = 100;  // t = 1.0

    for (int i = 0; i < n_steps; ++i) {
        euler.advance(U, dt, decay_rhs);
    }

    // Expected: exp(-1) ≈ 0.3679
    Real expected = std::exp(-1.0);
    // Euler is first order, expect O(dt) error
    EXPECT_NEAR(U[0].rho, expected, 0.02);
}

TEST_F(TimeIntegratorTest, SSPRK3ConvergesBetter) {
    ConservativeArray U(1, ConservativeVars{1.0, 0.0, 0.0});
    SSPRK3 rk3;

    Real dt = 0.01;
    int n_steps = 100;  // t = 1.0

    for (int i = 0; i < n_steps; ++i) {
        rk3.advance(U, dt, decay_rhs);
    }

    Real expected = std::exp(-1.0);
    // RK3 is third order, expect much smaller error
    EXPECT_NEAR(U[0].rho, expected, 1e-5);
}

TEST_F(TimeIntegratorTest, VariantDispatch) {
    ConservativeArray U(1, ConservativeVars{1.0, 0.0, 0.0});
    TimeIntegratorVariant integrator = SSPRK3{};

    advance(integrator, U, 0.01, decay_rhs);

    EXPECT_LT(U[0].rho, 1.0);  // Should have decayed
}

TEST_F(TimeIntegratorTest, MultipleVariables) {
    ConservativeArray U(10, ConservativeVars{1.0, 2.0, 3.0});
    SSPRK3 rk3;

    rk3.advance(U, 0.01, decay_rhs);

    for (const auto& u : U) {
        EXPECT_LT(u.rho, 1.0);
        EXPECT_LT(u.rho_u, 2.0);
        EXPECT_LT(u.E, 3.0);
    }
}
