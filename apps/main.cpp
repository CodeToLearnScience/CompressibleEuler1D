/**
 * @file main.cpp
 * @brief Entry point for the 1D Euler solver
 */

#include "euler1d/config/parser.hpp"
#include "euler1d/solver/solver.hpp"
#include "euler1d/io/output.hpp"

#include <filesystem>
#include <iostream>
#include <print>
#include <string>

void print_usage(const char* program) {
    std::println("Usage: {} <config.toml> [output_dir]", program);
    std::println("");
    std::println("Arguments:");
    std::println("  config.toml  Path to TOML configuration file");
    std::println("  output_dir   Optional output directory (default: current directory)");
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    const std::filesystem::path config_path{argv[1]};
    const std::filesystem::path output_dir = (argc >= 3) ? argv[2] : ".";

    // Create output directory if needed
    std::filesystem::create_directories(output_dir);

    try {
        // Parse configuration
        std::println("Loading configuration: {}", config_path.string());
        const auto config = euler1d::parse_config(config_path);

        // Create and run solver
        euler1d::Solver solver(config);
        solver.run();

        // Get solution
        const auto& U = solver.solution();
        const auto W = solver.to_primitive();
        const auto& mesh = solver.mesh();

        // Write output files
        const std::string base_name = config.simulation.test_name;

        const auto csv_path = output_dir / (base_name + ".csv");
        euler1d::write_csv(csv_path, mesh, U, W, solver.time());
        std::println("Wrote CSV: {}", csv_path.string());

        const auto vtk_path = output_dir / (base_name + ".vtk");
        euler1d::write_vtk(vtk_path, mesh, U, W, solver.time());
        std::println("Wrote VTK: {}", vtk_path.string());

        return 0;

    } catch (const euler1d::ConfigError& e) {
        std::println(stderr, "Configuration error: {}", e.what());
        return 1;
    } catch (const std::exception& e) {
        std::println(stderr, "Error: {}", e.what());
        return 1;
    }
}
