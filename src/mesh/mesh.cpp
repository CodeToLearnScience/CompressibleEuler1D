/**
 * @file mesh.cpp
 * @brief Mesh implementation
 */

#include "euler1d/mesh/mesh.hpp"
#include <stdexcept>

namespace euler1d {

Mesh1D::Mesh1D(Real xmin, Real xmax, int num_cells)
    : xmin_{xmin}, xmax_{xmax}, num_cells_{num_cells} {
    if (num_cells <= 0) {
        throw std::invalid_argument("num_cells must be positive");
    }
    if (xmax <= xmin) {
        throw std::invalid_argument("xmax must be greater than xmin");
    }
    dx_ = (xmax_ - xmin_) / static_cast<Real>(num_cells_);
}

}  // namespace euler1d
