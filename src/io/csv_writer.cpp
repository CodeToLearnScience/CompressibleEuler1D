/**
 * @file csv_writer.cpp
 * @brief CSV output implementation
 */

#include "euler1d/io/output.hpp"
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace euler1d {

void write_csv(const std::filesystem::path& path, const Mesh1D& mesh,
               const ConservativeArray& U, const PrimitiveArray& W, Real time) {
    std::ofstream file(path);
    if (!file) {
        throw std::runtime_error("Cannot open file for writing: " + path.string());
    }

    // Set precision
    file << std::setprecision(12) << std::scientific;

    // Header
    file << "# 1D Euler solution at time = " << time << "\n";
    file << "# x,rho,u,p,E\n";

    // Write interior cells only
    for (int i = mesh.first_interior(); i <= mesh.last_interior(); ++i) {
        const auto idx = static_cast<std::size_t>(i);
        file << mesh.x(i) << ","
             << W[idx].rho << ","
             << W[idx].u << ","
             << W[idx].p << ","
             << U[idx].E << "\n";
    }
}

}  // namespace euler1d
