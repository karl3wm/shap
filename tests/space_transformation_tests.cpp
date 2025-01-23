#include "shap/coord.hpp"
#include <shap/surface.hpp>
#include <shap/surfaces/flat_patch.hpp>
#include <shap/path.hpp>
#include "test_utils.hpp"
#include <cassert>
#include <iostream>

namespace shap::test {

// Test coordinate transformations with orthogonal basis
void test_orthogonal_basis_transformations() {
    auto face = surfaces::create_flat_patch(
        WorldPoint3(-1, 1, -1),
        WorldVector3(2, 0, 0),
        WorldVector3(0, 0, 2),
        1e-10,  // vector_length_epsilon
        1e-10   // parameter_bound_epsilon
    );

    // Test parameter to world mapping
    const std::vector<std::tuple<double, double, WorldPoint3>> param_points = {
        {0.0, 0.0, WorldPoint3(-1, 1, -1)},  // Origin
        {1.0, 0.0, WorldPoint3(1, 1, -1)},   // u = 1
        {0.0, 1.0, WorldPoint3(-1, 1, 1)},   // v = 1
        {0.5, 0.5, WorldPoint3(0, 1, 0)}     // Center
    };

    for (const auto& [u, v, expected] : param_points) {
        const auto geom = face->evaluate(ParamPoint2(u, v));
        assert(approx_equal(geom.world_pos(), expected));
    }

    // Test world to parameter mapping
    const std::vector<std::tuple<WorldPoint3, double, double>> world_points = {
        {WorldPoint3(-1, 1, -1), 0.0, 0.0},  // Origin
        {WorldPoint3(1, 1, -1), 1.0, 0.0},   // u = 1
        {WorldPoint3(-1, 1, 1), 0.0, 1.0},   // v = 1
        {WorldPoint3(0, 1, 0), 0.5, 0.5}     // Center
    };

    for (const auto& [world, expected_u, expected_v] : world_points) {
        const auto params = face->world_to_param(world).to_r2();
        assert(approx_equal(params.u(), expected_u) && 
               approx_equal(params.v(), expected_v));
    }
}

// Test points above and below surface
void test_normal_distance() {
    // Create a patch in the y=0 plane (normal along y axis)
    auto face = surfaces::create_flat_patch(
        WorldPoint3(0, 0, 0),      // origin at (0,0,0)
        WorldVector3(1, 0, 0),      // unit vector in x
        WorldVector3(0, 0, 1),      // unit vector in z
        1e-10,  // vector_length_epsilon
        1e-10   // parameter_bound_epsilon
    );

    // Point above surface (positive y)
    const WorldPoint3 above(0.5, 1.0, 0.5);  // Should map to u=0.5, v=0.5, normal=1.0
    const auto above_local = face->world_to_param(above);
    assert(approx_equal(above_local.u(), 0.5));
    assert(approx_equal(above_local.v(), 0.5));
    assert(approx_equal(above_local.w(), -1.0));  // Negative normal distance (normal points down)

    // Point below surface (negative y)
    const WorldPoint3 below(0.5, -1.0, 0.5);  // Should map to u=0.5, v=0.5, normal=-1.0
    const auto below_local = face->world_to_param(below);
    assert(approx_equal(below_local.u(), 0.5));
    assert(approx_equal(below_local.v(), 0.5));
    assert(approx_equal(below_local.w(), 1.0));  // Positive normal distance (normal points down)
}

// Test points outside parameter range
void test_outside_parameter_range() {
    auto face = surfaces::create_flat_patch(
        WorldPoint3(0, 0, 0),
        WorldVector3(1, 0, 0),
        WorldVector3(0, 0, 1),
        1e-10,
        1e-10
    );

    // Point beyond u=1 boundary
    const WorldPoint3 beyond_u(2.0, 0.0, 0.5);
    const auto beyond_u_local = face->world_to_param(beyond_u);
    assert(beyond_u_local.u() > 1.0);
    assert(approx_equal(beyond_u_local.v(), 0.5));

    // Point beyond v=1 boundary
    const WorldPoint3 beyond_v(0.5, 0.0, 2.0);
    const auto beyond_v_local = face->world_to_param(beyond_v);
    assert(approx_equal(beyond_v_local.u(), 0.5));
    assert(beyond_v_local.v() > 1.0);
}

// Test non-orthogonal basis
void test_skewed_basis() {
    auto face = surfaces::create_flat_patch(
        WorldPoint3(0, 0, 0),
        WorldVector3(1, 0, 0),          // First basis vector along x
        WorldVector3(0.5, 0, 1),        // Second basis vector skewed in x-z plane
        1e-10,
        1e-10
    );

    // Test point that should map to u=0.5, v=0.5
    const WorldPoint3 test_point(0.75, 0, 0.5);  // 0.75 = 0.5 + 0.5*0.5 (due to skew)
    const auto local = face->world_to_param(test_point);
    assert(approx_equal(local.u(), 0.5));
    assert(approx_equal(local.v(), 0.5));
}

// Test degenerate cases
void test_degenerate_cases() {
    try {
        // Create patch with nearly parallel basis vectors
        auto face = surfaces::create_flat_patch(
            WorldPoint3(0, 0, 0),
            WorldVector3(1, 0, 0),
            WorldVector3(1, 0, 1e-11),  // Almost parallel to first vector
            1e-10,
            1e-10
        );
        
        const WorldPoint3 test_point(0.5, 0, 0);
        [[maybe_unused]] const auto result = face->world_to_param(test_point);  // Should throw
        assert(false);  // Should not reach here
    }
    catch (const std::invalid_argument& e) {
        // Expected exception
        std::string error_msg = e.what();
        std::cout << "Caught expected error: " << error_msg << std::endl;
        assert(error_msg.find("parallel") != std::string::npos);
    }
}

// Test path creation and evaluation
void test_path_creation() {
    auto face = surfaces::create_flat_patch(
        WorldPoint3(-1, 1, -1),
        WorldVector3(2, 0, 0),
        WorldVector3(0, 0, 2),
        1e-10,  // vector_length_epsilon
        1e-10   // parameter_bound_epsilon
    );

    // Test path length preservation
    const WorldPoint3 start(-0.5, 1, 0);
    const WorldVector3 dir(1, 0, 0);
    const double length = 1.0;
    const WorldPoint3 expected_end(0.5, 1, 0);

    const auto params = face->world_to_param(start).to_r2();
    const auto start_point = face->evaluate(params);
    auto path = face->create_path(start_point, dir, length);

    // Check key points for distance preservation
    const std::vector<double> check_points = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    for (double t : check_points) {
        const auto pt = path->evaluate(t);
        const auto pos = pt.world_pos();
        const double actual_dist = (pos - start).length();
        const double expected_dist = t * length;

        std::cout << "t=" << t << " expected=" << expected_dist 
                 << " actual=" << actual_dist << "\n";
        assert(approx_equal(actual_dist, expected_dist));
    }

    // Verify end point
    const auto end_pt = path->evaluate(1.0);
    assert(approx_equal(end_pt.world_pos(), expected_end));
}

} // namespace shap::test

int main() {
    try {
        // Run all test cases
        shap::test::test_orthogonal_basis_transformations();
        shap::test::test_normal_distance();
        shap::test::test_outside_parameter_range();
        shap::test::test_skewed_basis();
        shap::test::test_degenerate_cases();
        shap::test::test_path_creation();
        
        std::cout << "All tests completed.\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << "\n";
        return 1;
    }
}
