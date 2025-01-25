#include "shap/coord.hpp"
#include <shap/surface3d.hpp>
#include <shap/surfaces/flat_patch.hpp>
#include <shap/geometric_point.hpp>
#include <shap/path.hpp>
#include "test_utils.hpp"
#include <cassert>
#include <iostream>

namespace shap::test {

// Tests that are currently passing
void test_path_length_invariants_passing() {
    auto face = surfaces::create_flat_patch(
        WorldPoint3(-1, 1, -1),
        WorldVector3(2, 0, 0),
        WorldVector3(0, 0, 2),
        1e-10,  // vector_length_epsilon
        1e-10   // parameter_bound_epsilon
    );

    // Test parameter space distance
    const WorldPoint3 start(-0.5, 1, 0);
    const WorldVector3 dir(1, 0, 0);
    const double length = 1.0;
    const double expected_param_delta = 0.25;  // L/(2|du|) = 1.0/(2*2)

    const auto params = face->nearest(start);
    const auto start_point = face->evaluate(params);
    auto path = face->create_path(start_point, dir, length);

    const auto end_pt = path->evaluate(ParamPoint1(1.0));
    const double actual_param_delta = end_pt.local_pos().u() - 
                                    start_point.local_pos().u();

    assert(approx_equal(actual_param_delta, expected_param_delta));
}

// Tests that are currently failing and under investigation
void test_path_length_invariants_failing() {
    auto face = surfaces::create_flat_patch(
        WorldPoint3(-1, 1, -1),
        WorldVector3(2, 0, 0),
        WorldVector3(0, 0, 2),
        1e-10,  // vector_length_epsilon
        1e-10   // parameter_bound_epsilon
    );

    // Test world space distance preservation
    const WorldPoint3 start(-0.5, 1, 0);
    const WorldVector3 dir(1, 0, 0);
    const double length = 1.0;

    const auto params = face->nearest(start);
    const auto start_point = face->evaluate(params);
    auto path = face->create_path(start_point, dir, length);

    // Check key points for distance preservation
    const std::vector<double> check_points = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    for (double t : check_points) {
        const auto pt = path->evaluate(ParamPoint1(t));
        const auto pos = pt.world_pos();
        const double actual_dist = (pos - start).length();
        const double expected_dist = t * length;

        std::cout << "t=" << t << " expected=" << expected_dist 
                 << " actual=" << actual_dist << "\n";
        assert(approx_equal(actual_dist, expected_dist));
    }
}

} // namespace shap::test

int main() {
    try {
        // Run passing tests first
        shap::test::test_path_length_invariants_passing();

        // Run failing tests separately
        std::cout << "\n----------------------------------------\n";
        std::cout << "Running tests with known failures:\n";
        std::cout << "----------------------------------------\n";
        shap::test::test_path_length_invariants_failing();
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << "\n";
        return 1;
    }
}
