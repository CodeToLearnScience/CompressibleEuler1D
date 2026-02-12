/**
 * @file test_initial_condition.cpp
 * @brief Unit tests for initial conditions
 */

#include <gtest/gtest.h>
#include "euler1d/initial/initial_condition.hpp"
#include "euler1d/eos/eos.hpp"
#include <cmath>

using namespace euler1d;

class InitialConditionTest : public ::testing::Test {
protected:
    Mesh1D mesh{0.0, 1.0, 100};
    IdealGas eos{1.4};
    ConservativeArray U;

    void SetUp() override {
        U.resize(static_cast<std::size_t>(mesh.total_cells()));
    }
};

TEST_F(InitialConditionTest, PiecewiseConstantTwoRegions) {
    PiecewiseConstantIC ic;
    ic.regions = {
        Region{0.0, 0.5, 1.0, 0.0, 1.0},    // Left half
        Region{0.5, 1.0, 0.125, 0.0, 0.1}   // Right half
    };

    ic.apply(U, mesh, eos);

    // Check a cell in left region
    int left_cell = mesh.first_interior() + 10;
    auto W_left = eos.to_primitive(U[static_cast<std::size_t>(left_cell)]);
    EXPECT_NEAR(W_left.rho, 1.0, 1e-10);
    EXPECT_NEAR(W_left.p, 1.0, 1e-10);

    // Check a cell in right region
    int right_cell = mesh.last_interior() - 10;
    auto W_right = eos.to_primitive(U[static_cast<std::size_t>(right_cell)]);
    EXPECT_NEAR(W_right.rho, 0.125, 1e-10);
    EXPECT_NEAR(W_right.p, 0.1, 1e-10);
}

TEST_F(InitialConditionTest, ShockEntropyInteraction) {
    Mesh1D mesh2{-1.0, 1.0, 200};
    U.resize(static_cast<std::size_t>(mesh2.total_cells()));

    ShockEntropyInteractionIC ic;
    ic.discontinuity_position = 0.0;
    ic.left_state = {1.0, 0.0, 1.0};
    ic.right_state = {1.0, 0.2, 5.0, true, 0.0, 1.0};  // rho = 1 + 0.2*sin(5*pi*x)

    ic.apply(U, mesh2, eos);

    // Left of discontinuity: constant
    int left_cell = mesh2.first_interior() + 10;  // x < 0
    auto W_left = eos.to_primitive(U[static_cast<std::size_t>(left_cell)]);
    EXPECT_NEAR(W_left.rho, 1.0, 1e-10);

    // Right of discontinuity: sinusoidal density
    int right_cell = mesh2.last_interior() - 10;  // x > 0
    auto W_right = eos.to_primitive(U[static_cast<std::size_t>(right_cell)]);
    Real x = mesh2.x(right_cell);
    Real expected_rho = 1.0 + 0.2 * std::sin(5.0 * constants::pi * x);
    EXPECT_NEAR(W_right.rho, expected_rho, 1e-10);
}

TEST_F(InitialConditionTest, CreateFromConfig) {
    InitialConditionConfig config;
    config.type = InitialConditionType::PiecewiseConstant;
    config.regions = {
        Region{0.0, 1.0, 1.0, 0.0, 1.0}
    };

    auto ic = create_initial_condition(config);

    // Should create PiecewiseConstantIC
    EXPECT_TRUE(std::holds_alternative<PiecewiseConstantIC>(ic));
}
