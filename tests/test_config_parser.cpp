/**
 * @file test_config_parser.cpp
 * @brief Unit tests for TOML configuration parser
 */

#include <gtest/gtest.h>
#include "euler1d/config/parser.hpp"
#include <filesystem>

using namespace euler1d;

class ConfigParserTest : public ::testing::Test {
protected:
    std::filesystem::path data_dir{"data"};
};

TEST_F(ConfigParserTest, ParseTestCase1) {
    auto config = parse_config(data_dir / "test_case1.toml");

    EXPECT_EQ(config.simulation.test_name, "test_case1");
    EXPECT_EQ(config.mesh.num_cells, 1000);
    EXPECT_DOUBLE_EQ(config.mesh.xmin, 0.0);
    EXPECT_DOUBLE_EQ(config.mesh.xmax, 1.0);
    EXPECT_DOUBLE_EQ(config.time.cfl, 0.5);
    EXPECT_DOUBLE_EQ(config.time.final_time, 0.2);
    EXPECT_EQ(config.time.integrator, TimeIntegrator::SSPRK3);
    EXPECT_EQ(config.numerics.order, 1);
    EXPECT_EQ(config.numerics.flux, FluxScheme::LLF);
    EXPECT_EQ(config.numerics.limiter, Limiter::VanLeer);
    EXPECT_DOUBLE_EQ(config.eos.gamma, 1.4);
    EXPECT_EQ(config.boundary.left, BoundaryType::Transmissive);
    EXPECT_EQ(config.boundary.right, BoundaryType::Transmissive);

    // Initial condition
    EXPECT_EQ(config.initial_condition.type, InitialConditionType::PiecewiseConstant);
    ASSERT_EQ(config.initial_condition.regions.size(), 2);

    // First region
    EXPECT_DOUBLE_EQ(config.initial_condition.regions[0].x_left, 0.0);
    EXPECT_DOUBLE_EQ(config.initial_condition.regions[0].x_right, 0.3);
    EXPECT_DOUBLE_EQ(config.initial_condition.regions[0].rho, 1.0);
    EXPECT_DOUBLE_EQ(config.initial_condition.regions[0].u, 0.75);
    EXPECT_DOUBLE_EQ(config.initial_condition.regions[0].p, 1.0);

    // Second region
    EXPECT_DOUBLE_EQ(config.initial_condition.regions[1].x_left, 0.3);
    EXPECT_DOUBLE_EQ(config.initial_condition.regions[1].x_right, 1.0);
    EXPECT_DOUBLE_EQ(config.initial_condition.regions[1].rho, 0.125);
    EXPECT_DOUBLE_EQ(config.initial_condition.regions[1].u, 0.0);
    EXPECT_DOUBLE_EQ(config.initial_condition.regions[1].p, 0.1);
}

TEST_F(ConfigParserTest, ParseTestCase10WithThreeRegions) {
    auto config = parse_config(data_dir / "test_case10.toml");

    EXPECT_EQ(config.initial_condition.type, InitialConditionType::PiecewiseConstant);
    ASSERT_EQ(config.initial_condition.regions.size(), 3);

    // Check all three regions exist
    EXPECT_DOUBLE_EQ(config.initial_condition.regions[0].p, 1000.0);
    EXPECT_DOUBLE_EQ(config.initial_condition.regions[1].p, 0.01);
    EXPECT_DOUBLE_EQ(config.initial_condition.regions[2].p, 100.0);
}

TEST_F(ConfigParserTest, ParseTestCase11ShockEntropy) {
    auto config = parse_config(data_dir / "test_case11.toml");

    EXPECT_EQ(config.initial_condition.type, InitialConditionType::ShockEntropyInteraction);
    EXPECT_DOUBLE_EQ(config.initial_condition.discontinuity_position, -0.8);

    // Left state
    EXPECT_DOUBLE_EQ(config.initial_condition.left_state.rho, 3.857143);
    EXPECT_DOUBLE_EQ(config.initial_condition.left_state.u, 2.629369);
    EXPECT_DOUBLE_EQ(config.initial_condition.left_state.p, 10.33333);

    // Right state (sinusoidal)
    EXPECT_DOUBLE_EQ(config.initial_condition.right_state.rho_base, 1.0);
    EXPECT_DOUBLE_EQ(config.initial_condition.right_state.rho_amplitude, 0.2);
    EXPECT_DOUBLE_EQ(config.initial_condition.right_state.rho_frequency, 5.0);
    EXPECT_TRUE(config.initial_condition.right_state.use_pi);  // "pi" function
}

TEST_F(ConfigParserTest, AllTestCasesParseSuccessfully) {
    for (int i = 1; i <= 12; ++i) {
        std::string filename = "test_case" + std::to_string(i) + ".toml";
        EXPECT_NO_THROW(parse_config(data_dir / filename)) << "Failed to parse " << filename;
    }
}

TEST_F(ConfigParserTest, InvalidFileThrows) {
    EXPECT_THROW(parse_config("nonexistent.toml"), ConfigError);
}
