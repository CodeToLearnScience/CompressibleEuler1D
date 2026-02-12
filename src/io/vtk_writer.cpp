/**
 * @file vtk_writer.cpp
 * @brief VTK output implementation
 */

#include "euler1d/io/output.hpp"
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace euler1d {

void write_vtk(const std::filesystem::path& path, const Mesh1D& mesh,
               const ConservativeArray& U, const PrimitiveArray& W, Real /*time*/) {
    std::ofstream file(path);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + path.string());
    }

    const int n = mesh.num_cells();
    const int first = mesh.first_interior();

    file << std::setprecision(12) << std::scientific;

    // VTK header
    file << "# vtk DataFile Version 3.0\n";
    file << "1D Euler solution\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << n << " 1 1\n";
    file << "POINTS " << n << " double\n";

    // Point coordinates (cell centers)
    for (int i = 0; i < n; ++i) {
        file << mesh.x(first + i) << " 0 0\n";
    }

    // Cell data
    file << "\nPOINT_DATA " << n << "\n";

    // Density
    file << "SCALARS rho double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < n; ++i) {
        file << W[static_cast<std::size_t>(first + i)].rho << "\n";
    }

    // Velocity
    file << "\nSCALARS u double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < n; ++i) {
        file << W[static_cast<std::size_t>(first + i)].u << "\n";
    }

    // Pressure
    file << "\nSCALARS p double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < n; ++i) {
        file << W[static_cast<std::size_t>(first + i)].p << "\n";
    }

    // Total energy
    file << "\nSCALARS E double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < n; ++i) {
        file << U[static_cast<std::size_t>(first + i)].E << "\n";
    }
}

}  // namespace euler1d
