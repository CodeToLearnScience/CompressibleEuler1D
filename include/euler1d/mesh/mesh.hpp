/**
 * @file mesh.hpp
 * @brief 1D mesh with ghost cells for the Euler solver
 */

#ifndef EULER1D_MESH_MESH_HPP
#define EULER1D_MESH_MESH_HPP

#include "../core/types.hpp"
#include <cstddef>

namespace euler1d {

/**
 * @brief 1D uniform mesh with ghost cells
 *
 * The mesh is structured as:
 *   [ghost][ghost] | [interior cells] | [ghost][ghost]
 *   0      1        2 ... num_cells+1   num_cells+2  num_cells+3
 *
 * For second-order schemes, we use 2 ghost cells on each side.
 */
class Mesh1D {
public:
    static constexpr int num_ghosts = 2;  ///< Ghost cells per side

    /// Construct mesh from physical domain and cell count
    Mesh1D(Real xmin, Real xmax, int num_cells);

    /// Physical domain bounds
    [[nodiscard]] Real xmin() const noexcept { return xmin_; }
    [[nodiscard]] Real xmax() const noexcept { return xmax_; }

    /// Number of interior cells (excluding ghosts)
    [[nodiscard]] int num_cells() const noexcept { return num_cells_; }

    /// Total number of cells (including ghosts)
    [[nodiscard]] int total_cells() const noexcept { return num_cells_ + 2 * num_ghosts; }

    /// Cell size (uniform)
    [[nodiscard]] Real dx() const noexcept { return dx_; }

    /// Cell center coordinate for cell i (including ghost cells)
    [[nodiscard]] Real x(int i) const noexcept {
        // Ghost cells extend beyond domain
        return xmin_ + (static_cast<Real>(i) - num_ghosts + Real{0.5}) * dx_;
    }

    /// Left face coordinate for cell i
    [[nodiscard]] Real x_face_left(int i) const noexcept {
        return xmin_ + (static_cast<Real>(i) - num_ghosts) * dx_;
    }

    /// Right face coordinate for cell i
    [[nodiscard]] Real x_face_right(int i) const noexcept {
        return x_face_left(i) + dx_;
    }

    /// First interior cell index
    [[nodiscard]] static constexpr int first_interior() noexcept { return num_ghosts; }

    /// Last interior cell index (inclusive)
    [[nodiscard]] int last_interior() const noexcept { return num_ghosts + num_cells_ - 1; }

    /// Check if index is an interior cell
    [[nodiscard]] bool is_interior(int i) const noexcept {
        return i >= first_interior() && i <= last_interior();
    }

private:
    Real xmin_;
    Real xmax_;
    int num_cells_;
    Real dx_;
};

}  // namespace euler1d

#endif  // EULER1D_MESH_MESH_HPP
