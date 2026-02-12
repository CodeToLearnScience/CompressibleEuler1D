/**
 * @file parser.cpp
 * @brief TOML configuration parser implementation
 */

#include "euler1d/config/parser.hpp"
#include <toml++/toml.hpp>
#include <algorithm>
#include <cctype>

namespace euler1d {

namespace {

/// Convert string to lowercase
std::string to_lower(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return str;
}

}  // namespace

// =============================================================================
// String to enum conversions
// =============================================================================

FluxScheme parse_flux_scheme(const std::string& str) {
    const auto lower = to_lower(str);
    if (lower == "llf" || lower == "local_lax_friedrichs") return FluxScheme::LLF;
    if (lower == "rusanov") return FluxScheme::Rusanov;
    if (lower == "hll") return FluxScheme::HLL;
    if (lower == "hllc") return FluxScheme::HLLC;
    throw ConfigError("Unknown flux scheme: " + str);
}

Limiter parse_limiter(const std::string& str) {
    const auto lower = to_lower(str);
    if (lower == "none" || lower == "nolimiter") return Limiter::None;
    if (lower == "minmod") return Limiter::Minmod;
    if (lower == "vanleer" || lower == "van_leer") return Limiter::VanLeer;
    if (lower == "superbee") return Limiter::Superbee;
    if (lower == "mc" || lower == "monotonized_central") return Limiter::MC;
    throw ConfigError("Unknown limiter: " + str);
}

TimeIntegrator parse_time_integrator(const std::string& str) {
    const auto lower = to_lower(str);
    if (lower == "euler" || lower == "explicit_euler" || lower == "forward_euler") {
        return TimeIntegrator::ExplicitEuler;
    }
    if (lower == "ssprk3" || lower == "rk3" || lower == "ssp_rk3") {
        return TimeIntegrator::SSPRK3;
    }
    throw ConfigError("Unknown time integrator: " + str);
}

BoundaryType parse_boundary_type(const std::string& str) {
    const auto lower = to_lower(str);
    if (lower == "transmissive" || lower == "outflow" || lower == "zero_gradient") {
        return BoundaryType::Transmissive;
    }
    if (lower == "reflective" || lower == "wall" || lower == "solid_wall") {
        return BoundaryType::Reflective;
    }
    if (lower == "periodic") return BoundaryType::Periodic;
    throw ConfigError("Unknown boundary type: " + str);
}

EosModel parse_eos_model(const std::string& str) {
    const auto lower = to_lower(str);
    if (lower == "ideal_gas" || lower == "idealgas") return EosModel::IdealGas;
    throw ConfigError("Unknown EOS model: " + str);
}

InitialConditionType parse_initial_condition_type(const std::string& str) {
    const auto lower = to_lower(str);
    if (lower == "piecewise_constant" || lower == "piecewiseconstant") {
        return InitialConditionType::PiecewiseConstant;
    }
    if (lower == "shock_entropy_interaction" || lower == "shockentropyinteraction" ||
        lower == "shock_entropy" || lower == "shu_osher") {
        return InitialConditionType::ShockEntropyInteraction;
    }
    throw ConfigError("Unknown initial condition type: " + str);
}

// =============================================================================
// Main parser
// =============================================================================

Config parse_config(const std::filesystem::path& path) {
    Config config;

    toml::table tbl;
    try {
        tbl = toml::parse_file(path.string());
    } catch (const toml::parse_error& e) {
        throw ConfigError("TOML parse error: " + std::string(e.description()));
    }

    // [simulation]
    if (auto sim = tbl["simulation"].as_table()) {
        if (auto v = (*sim)["equations"].value<std::string>()) {
            config.simulation.equations = *v;
        }
        if (auto v = (*sim)["test_name"].value<std::string>()) {
            config.simulation.test_name = *v;
        }
    }

    // [mesh]
    if (auto mesh = tbl["mesh"].as_table()) {
        if (auto v = (*mesh)["xmin"].value<double>()) {
            config.mesh.xmin = static_cast<Real>(*v);
        }
        if (auto v = (*mesh)["xmax"].value<double>()) {
            config.mesh.xmax = static_cast<Real>(*v);
        }
        if (auto v = (*mesh)["num_cells"].value<int64_t>()) {
            config.mesh.num_cells = static_cast<int>(*v);
        }
    }

    // [time]
    if (auto time = tbl["time"].as_table()) {
        if (auto v = (*time)["cfl"].value<double>()) {
            config.time.cfl = static_cast<Real>(*v);
        }
        if (auto v = (*time)["final_time"].value<double>()) {
            config.time.final_time = static_cast<Real>(*v);
        }
        if (auto v = (*time)["time_integrator"].value<std::string>()) {
            config.time.integrator = parse_time_integrator(*v);
        }
    }

    // [numerics]
    if (auto num = tbl["numerics"].as_table()) {
        if (auto v = (*num)["order"].value<int64_t>()) {
            config.numerics.order = static_cast<int>(*v);
        }
        if (auto v = (*num)["flux"].value<std::string>()) {
            config.numerics.flux = parse_flux_scheme(*v);
        }
        if (auto v = (*num)["limiter"].value<std::string>()) {
            config.numerics.limiter = parse_limiter(*v);
        }
    }

    // [eos]
    if (auto eos = tbl["eos"].as_table()) {
        if (auto v = (*eos)["model"].value<std::string>()) {
            config.eos.model = parse_eos_model(*v);
        }
        if (auto v = (*eos)["gamma"].value<double>()) {
            config.eos.gamma = static_cast<Real>(*v);
        }
    }

    // [boundary_conditions]
    if (auto bc = tbl["boundary_conditions"].as_table()) {
        if (auto v = (*bc)["left"].value<std::string>()) {
            config.boundary.left = parse_boundary_type(*v);
        }
        if (auto v = (*bc)["right"].value<std::string>()) {
            config.boundary.right = parse_boundary_type(*v);
        }
    }

    // [initial_condition]
    if (auto ic = tbl["initial_condition"].as_table()) {
        if (auto v = (*ic)["type"].value<std::string>()) {
            config.initial_condition.type = parse_initial_condition_type(*v);
        }

        // Parse based on type
        if (config.initial_condition.type == InitialConditionType::PiecewiseConstant) {
            // Parse [[initial_condition.region]] array
            if (auto regions = (*ic)["region"].as_array()) {
                for (const auto& elem : *regions) {
                    if (auto reg = elem.as_table()) {
                        Region region;
                        if (auto v = (*reg)["x_left"].value<double>()) {
                            region.x_left = static_cast<Real>(*v);
                        }
                        if (auto v = (*reg)["x_right"].value<double>()) {
                            region.x_right = static_cast<Real>(*v);
                        }
                        if (auto v = (*reg)["rho"].value<double>()) {
                            region.rho = static_cast<Real>(*v);
                        }
                        if (auto v = (*reg)["u"].value<double>()) {
                            region.u = static_cast<Real>(*v);
                        }
                        if (auto v = (*reg)["p"].value<double>()) {
                            region.p = static_cast<Real>(*v);
                        }
                        config.initial_condition.regions.push_back(region);
                    }
                }
            }
        } else if (config.initial_condition.type == InitialConditionType::ShockEntropyInteraction) {
            // Parse discontinuity position
            if (auto v = (*ic)["discontinuity_position"].value<double>()) {
                config.initial_condition.discontinuity_position = static_cast<Real>(*v);
            }

            // Parse left state (constant)
            if (auto left = (*ic)["left_state"].as_table()) {
                if (auto v = (*left)["rho"].value<double>()) {
                    config.initial_condition.left_state.rho = static_cast<Real>(*v);
                }
                if (auto v = (*left)["u"].value<double>()) {
                    config.initial_condition.left_state.u = static_cast<Real>(*v);
                }
                if (auto v = (*left)["p"].value<double>()) {
                    config.initial_condition.left_state.p = static_cast<Real>(*v);
                }
            }

            // Parse right state (sinusoidal)
            if (auto right = (*ic)["right_state"].as_table()) {
                if (auto v = (*right)["rho_base"].value<double>()) {
                    config.initial_condition.right_state.rho_base = static_cast<Real>(*v);
                }
                if (auto v = (*right)["rho_amplitude"].value<double>()) {
                    config.initial_condition.right_state.rho_amplitude = static_cast<Real>(*v);
                }
                if (auto v = (*right)["rho_frequency"].value<double>()) {
                    config.initial_condition.right_state.rho_frequency = static_cast<Real>(*v);
                }
                if (auto v = (*right)["rho_function"].value<std::string>()) {
                    const auto lower = to_lower(*v);
                    config.initial_condition.right_state.use_pi = (lower == "pi");
                }
                if (auto v = (*right)["u"].value<double>()) {
                    config.initial_condition.right_state.u = static_cast<Real>(*v);
                }
                if (auto v = (*right)["p"].value<double>()) {
                    config.initial_condition.right_state.p = static_cast<Real>(*v);
                }
            }
        }
    }

    return config;
}

}  // namespace euler1d
