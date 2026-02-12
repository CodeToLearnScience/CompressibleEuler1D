/**
 * @file boundary.hpp
 * @brief Boundary conditions for the 1D Euler solver
 */

#ifndef EULER1D_BOUNDARY_BOUNDARY_HPP
#define EULER1D_BOUNDARY_BOUNDARY_HPP

#include "../core/types.hpp"
#include "../mesh/mesh.hpp"
#include <span>
#include <variant>

namespace euler1d {

// =============================================================================
// Transmissive (Zero-Gradient) Boundary
// =============================================================================

/**
 * @brief Transmissive boundary condition (outflow/zero-gradient)
 *
 * Copies interior values to ghost cells. Allows waves to exit cleanly.
 */
struct TransmissiveBoundary {
    void apply_left(std::span<ConservativeVars> U, const Mesh1D& mesh) const {
        const int first = mesh.first_interior();
        for (int i = first - 1; i >= 0; --i) {
            U[static_cast<std::size_t>(i)] = U[static_cast<std::size_t>(first)];
        }
    }

    void apply_right(std::span<ConservativeVars> U, const Mesh1D& mesh) const {
        const int last = mesh.last_interior();
        const int n = mesh.total_cells();
        for (int i = last + 1; i < n; ++i) {
            U[static_cast<std::size_t>(i)] = U[static_cast<std::size_t>(last)];
        }
    }
};

// =============================================================================
// Reflective (Wall) Boundary
// =============================================================================

/**
 * @brief Reflective boundary condition (solid wall)
 *
 * Reflects the velocity component normal to the wall.
 */
struct ReflectiveBoundary {
    void apply_left(std::span<ConservativeVars> U, const Mesh1D& mesh) const {
        const int first = mesh.first_interior();
        for (int g = 0; g < Mesh1D::num_ghosts; ++g) {
            const int ghost_idx = first - 1 - g;
            const int interior_idx = first + g;
            auto& ghost = U[static_cast<std::size_t>(ghost_idx)];
            const auto& interior = U[static_cast<std::size_t>(interior_idx)];

            ghost.rho = interior.rho;
            ghost.rho_u = -interior.rho_u;  // Reflect velocity
            ghost.E = interior.E;
        }
    }

    void apply_right(std::span<ConservativeVars> U, const Mesh1D& mesh) const {
        const int last = mesh.last_interior();
        for (int g = 0; g < Mesh1D::num_ghosts; ++g) {
            const int ghost_idx = last + 1 + g;
            const int interior_idx = last - g;
            auto& ghost = U[static_cast<std::size_t>(ghost_idx)];
            const auto& interior = U[static_cast<std::size_t>(interior_idx)];

            ghost.rho = interior.rho;
            ghost.rho_u = -interior.rho_u;  // Reflect velocity
            ghost.E = interior.E;
        }
    }
};

// =============================================================================
// Periodic Boundary
// =============================================================================

/**
 * @brief Periodic boundary condition
 *
 * Left ghosts = right interior, right ghosts = left interior.
 */
struct PeriodicBoundary {
    void apply_left(std::span<ConservativeVars> U, const Mesh1D& mesh) const {
        const int first = mesh.first_interior();
        const int last = mesh.last_interior();
        for (int g = 0; g < Mesh1D::num_ghosts; ++g) {
            const int ghost_idx = first - 1 - g;
            const int source_idx = last - g;
            U[static_cast<std::size_t>(ghost_idx)] = U[static_cast<std::size_t>(source_idx)];
        }
    }

    void apply_right(std::span<ConservativeVars> U, const Mesh1D& mesh) const {
        const int first = mesh.first_interior();
        const int last = mesh.last_interior();
        for (int g = 0; g < Mesh1D::num_ghosts; ++g) {
            const int ghost_idx = last + 1 + g;
            const int source_idx = first + g;
            U[static_cast<std::size_t>(ghost_idx)] = U[static_cast<std::size_t>(source_idx)];
        }
    }
};

// =============================================================================
// Boundary Variant for runtime selection
// =============================================================================

/// Variant holding all supported boundary conditions
using BoundaryVariant = std::variant<TransmissiveBoundary, ReflectiveBoundary, PeriodicBoundary>;

/// Apply left boundary condition
inline void apply_left_boundary(const BoundaryVariant& bc, std::span<ConservativeVars> U, const Mesh1D& mesh) {
    std::visit([&](const auto& b) { b.apply_left(U, mesh); }, bc);
}

/// Apply right boundary condition
inline void apply_right_boundary(const BoundaryVariant& bc, std::span<ConservativeVars> U, const Mesh1D& mesh) {
    std::visit([&](const auto& b) { b.apply_right(U, mesh); }, bc);
}

/// Apply both boundary conditions
inline void apply_boundaries(const BoundaryVariant& bc_left, const BoundaryVariant& bc_right,
                              std::span<ConservativeVars> U, const Mesh1D& mesh) {
    apply_left_boundary(bc_left, U, mesh);
    apply_right_boundary(bc_right, U, mesh);
}

}  // namespace euler1d

#endif  // EULER1D_BOUNDARY_BOUNDARY_HPP
