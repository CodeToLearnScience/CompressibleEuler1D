/**
 * @file test_solver_integration.cpp
 * @brief Integration tests for the complete solver
 */

#include <gtest/gtest.h>
#include "euler1d/config/parser.hpp"
#include "euler1d/solver/solver.hpp"
#include <filesystem>
#include <cmath>

using namespace euler1d;

class SolverIntegrationTest : public ::testing::Test {
protected:
    std::filesystem::path data_dir{"data"};
};

TEST_F(SolverIntegrationTest, SodShockTubeRuns) {
    // Use test_case1 (Sod-like shock tube)
    auto config = parse_config(data_dir / "test_case1.toml");

    // Reduce final time for faster test
    config.time.final_time = 0.01;

    Solver solver(config);
    EXPECT_NO_THROW(solver.run());

    // Check solution is valid
    const auto& U = solver.solution();
    for (const auto& u : U) {
        EXPECT_TRUE(std::isfinite(u.rho));
        EXPECT_TRUE(std::isfinite(u.rho_u));
        EXPECT_TRUE(std::isfinite(u.E));
        EXPECT_GT(u.rho, 0.0);  // Positive density
    }
}

TEST_F(SolverIntegrationTest, ConservationOfMass) {
    auto config = parse_config(data_dir / "test_case1.toml");
    config.mesh.num_cells = 100;  // Smaller for faster test
    config.time.final_time = 0.05;

    Solver solver(config);

    // Compute initial mass
    Real initial_mass = 0;
    for (int i = solver.mesh().first_interior(); i <= solver.mesh().last_interior(); ++i) {
        initial_mass += solver.solution()[static_cast<std::size_t>(i)].rho;
    }
    initial_mass *= solver.mesh().dx();

    solver.run();

    // Compute final mass
    Real final_mass = 0;
    for (int i = solver.mesh().first_interior(); i <= solver.mesh().last_interior(); ++i) {
        final_mass += solver.solution()[static_cast<std::size_t>(i)].rho;
    }
    final_mass *= solver.mesh().dx();

    // Mass should be conserved (with some tolerance due to boundary conditions)
    EXPECT_NEAR(final_mass, initial_mass, 0.1);
}

TEST_F(SolverIntegrationTest, PositivityPreserved) {
    // Strong shock case
    auto config = parse_config(data_dir / "test_case3.toml");
    config.mesh.num_cells = 200;
    config.time.final_time = 0.005;

    Solver solver(config);
    solver.run();

    const auto W = solver.to_primitive();
    for (int i = solver.mesh().first_interior(); i <= solver.mesh().last_interior(); ++i) {
        EXPECT_GT(W[static_cast<std::size_t>(i)].rho, 0.0) << "Negative density at cell " << i;
        EXPECT_GT(W[static_cast<std::size_t>(i)].p, 0.0) << "Negative pressure at cell " << i;
    }
}

TEST_F(SolverIntegrationTest, ShockEntropyWaveRuns) {
    // Test the sinusoidal initial condition
    auto config = parse_config(data_dir / "test_case11.toml");
    config.mesh.num_cells = 100;
    config.time.final_time = 0.01;

    Solver solver(config);
    EXPECT_NO_THROW(solver.run());
}
