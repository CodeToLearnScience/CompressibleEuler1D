/**
 * @file output.hpp
 * @brief Output writers for the 1D Euler solver
 */

#ifndef EULER1D_IO_OUTPUT_HPP
#define EULER1D_IO_OUTPUT_HPP

#include "../core/types.hpp"
#include "../mesh/mesh.hpp"
#include <filesystem>
#include <string>

namespace euler1d {

/**
 * @brief Write solution to CSV file
 *
 * Columns: x, rho, u, p, E
 *
 * @param path Output file path
 * @param mesh Computational mesh
 * @param U Conservative solution
 * @param W Primitive solution
 * @param time Current simulation time
 */
void write_csv(const std::filesystem::path& path, const Mesh1D& mesh,
               const ConservativeArray& U, const PrimitiveArray& W, Real time);

/**
 * @brief Write solution to VTK legacy format
 *
 * Compatible with ParaView for visualization.
 *
 * @param path Output file path
 * @param mesh Computational mesh
 * @param U Conservative solution
 * @param W Primitive solution
 * @param time Current simulation time
 */
void write_vtk(const std::filesystem::path& path, const Mesh1D& mesh,
               const ConservativeArray& U, const PrimitiveArray& W, Real time);

}  // namespace euler1d

#endif  // EULER1D_IO_OUTPUT_HPP
