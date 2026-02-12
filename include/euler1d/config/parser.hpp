/**
 * @file parser.hpp
 * @brief TOML configuration parser
 */

#ifndef EULER1D_CONFIG_PARSER_HPP
#define EULER1D_CONFIG_PARSER_HPP

#include "config_types.hpp"
#include <filesystem>
#include <stdexcept>

namespace euler1d {

/// Exception thrown for configuration parsing errors
class ConfigError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/**
 * @brief Parse a TOML configuration file
 *
 * @param path Path to the TOML file
 * @return Parsed configuration
 * @throws ConfigError if parsing fails or required fields are missing
 */
Config parse_config(const std::filesystem::path& path);

}  // namespace euler1d

#endif  // EULER1D_CONFIG_PARSER_HPP
