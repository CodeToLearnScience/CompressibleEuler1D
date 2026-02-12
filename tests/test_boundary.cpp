/**
 * @file test_boundary.cpp
 * @brief Unit tests for boundary conditions
 */

#include <gtest/gtest.h>
#include "euler1d/boundary/boundary.hpp"

using namespace euler1d;

class BoundaryTest : public ::testing::Test {
protected:
    Mesh1D mesh{0.0, 1.0, 10};  // 10 interior cells + 4 ghost cells
    ConservativeArray U;

    void SetUp() override {
        U.resize(static_cast<std::size_t>(mesh.total_cells()));
        // Initialize interior cells
        for (int i = mesh.first_interior(); i <= mesh.last_interior(); ++i) {
            U[static_cast<std::size_t>(i)] = ConservativeVars{1.0, 2.0, 3.0};
        }
        // Set different values for first and last interior
        U[static_cast<std::size_t>(mesh.first_interior())] = ConservativeVars{1.0, 1.0, 1.0};
        U[static_cast<std::size_t>(mesh.last_interior())] = ConservativeVars{2.0, 2.0, 2.0};
    }
};

TEST_F(BoundaryTest, TransmissiveLeft) {
    TransmissiveBoundary bc;
    bc.apply_left(U, mesh);

    // Ghost cells should equal first interior cell
    const auto& first = U[static_cast<std::size_t>(mesh.first_interior())];
    EXPECT_DOUBLE_EQ(U[0].rho, first.rho);
    EXPECT_DOUBLE_EQ(U[1].rho, first.rho);
}

TEST_F(BoundaryTest, TransmissiveRight) {
    TransmissiveBoundary bc;
    bc.apply_right(U, mesh);

    // Ghost cells should equal last interior cell
    const auto& last = U[static_cast<std::size_t>(mesh.last_interior())];
    const int n = mesh.total_cells();
    EXPECT_DOUBLE_EQ(U[static_cast<std::size_t>(n - 2)].rho, last.rho);
    EXPECT_DOUBLE_EQ(U[static_cast<std::size_t>(n - 1)].rho, last.rho);
}

TEST_F(BoundaryTest, ReflectiveLeft) {
    ReflectiveBoundary bc;
    bc.apply_left(U, mesh);

    // Velocity should be reflected
    const auto& first = U[static_cast<std::size_t>(mesh.first_interior())];
    EXPECT_DOUBLE_EQ(U[1].rho, first.rho);
    EXPECT_DOUBLE_EQ(U[1].rho_u, -first.rho_u);  // Reflected
    EXPECT_DOUBLE_EQ(U[1].E, first.E);
}

TEST_F(BoundaryTest, PeriodicBoundary) {
    PeriodicBoundary bc;

    // Apply both boundaries
    bc.apply_left(U, mesh);
    bc.apply_right(U, mesh);

    // Left ghosts should match right interior
    const auto& last = U[static_cast<std::size_t>(mesh.last_interior())];
    EXPECT_DOUBLE_EQ(U[1].rho, last.rho);

    // Right ghosts should match left interior
    const auto& first = U[static_cast<std::size_t>(mesh.first_interior())];
    int n = mesh.total_cells();
    EXPECT_DOUBLE_EQ(U[static_cast<std::size_t>(n - 2)].rho, first.rho);
}

TEST_F(BoundaryTest, VariantDispatch) {
    BoundaryVariant bc = TransmissiveBoundary{};

    apply_left_boundary(bc, U, mesh);
    apply_right_boundary(bc, U, mesh);

    // Should not throw
    SUCCEED();
}
