#include "shap/coord.hpp"
#include <shap/surface.hpp>
#include <shap/surfaces/basic_surfaces.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <functional>

namespace shap::test {

// Test utilities
namespace {
    constexpr double EPSILON = 1e-10;
    
    // Helper to check if two doubles are approximately equal
    [[nodiscard]] constexpr bool approx_equal(double a, double b, double epsilon = EPSILON) noexcept {
        return std::abs(a - b) <= epsilon;
    }

    // Helper to check if two points are approximately equal
    [[nodiscard]] bool approx_equal(const Point& a, const Point& b, double epsilon = EPSILON) noexcept {
        return approx_equal(a.x, b.x, epsilon) &&
               approx_equal(a.y, b.y, epsilon) &&
               approx_equal(a.z, b.z, epsilon);
    }

    // Print a point for debugging
    void print_point(std::string_view label, const Point& p) {
        std::cout << label << ": ("
                 << std::fixed << std::setprecision(6)
                 << p.x << ", " << p.y << ", " << p.z << ")\n";
    }

    // Test fixture for cube face tests
    class CubeFaceTest {
    public:
        CubeFaceTest() : face_(create_front_face()) {}

        // Create front face of unit cube (y = 1)
        static std::shared_ptr<surfaces::SquareFace> create_front_face() {
            return std::make_shared<surfaces::SquareFace>(
                Point(-1, 1, -1),    // origin at top-left
                Vector(2, 0, 0),     // u: left to right
                Vector(0, 0, 2)      // v: top to bottom
            );
        }

        // Helper to convert world coordinates to parameters
        [[nodiscard]] static ParameterCoordinate world_to_params(const Point& world) noexcept {
            return ParameterCoordinate{
                (world.x + 1) / 2,  // u = (x + 1)/2
                (world.z + 1) / 2   // v = (z + 1)/2
            };
        }

        std::shared_ptr<surfaces::SquareFace> face_;
    };

    // Test fixture for cube tests
    class CubeTest {
    public:
        CubeTest() : cube_(surfaces::create_cube(1.0)) {}
        SurfaceCollection cube_;
    };

    // Timing utility
    template<typename F>
    double time_operation(F&& func) {
        const auto start = std::chrono::high_resolution_clock::now();
        func();
        const auto end = std::chrono::high_resolution_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        return static_cast<double>(duration.count()) / 1000.0;  // Convert to milliseconds
    }
}

// Test straight line paths on cube face
void test_cube_face_paths() {
    CubeFaceTest test;
    
    // Test parallel path
    {
        const Point start(-0.5, 1, 0);
        const Vector dir(1, 0, 0);
        const double length = 1.0;
        const Point expected_end(0.5, 1, 0);

        const auto params = test.world_to_params(start);
        const auto start_point = test.face_->evaluate(params);
        auto path = test.face_->create_path(start_point, dir, length);

        // Check key points
        const std::vector<double> check_points = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
        for (double t : check_points) {
            const auto pt = path->evaluate(t);
            const auto pos = pt.world_pos();
            std::cout << "t=" << t << " pos=(" << pos.x << "," 
                     << pos.y << "," << pos.z << ")\n";
        }

        const auto end_pt = path->evaluate(1.0);
        assert(approx_equal(end_pt.world_pos(), expected_end));
    }
    
    // Test diagonal path
    {
        const Point start(-0.5, 1, -0.5);
        const Vector dir = Vector(1, 0, 1).normalized();
        const double length = std::sqrt(2);
        const Point expected_end(0.5, 1, 0.5);

        const auto params = test.world_to_params(start);
        const auto start_point = test.face_->evaluate(params);
        auto path = test.face_->create_path(start_point, dir, length);

        const auto end_pt = path->evaluate(length);
        assert(approx_equal(end_pt.world_pos(), expected_end));
    }

    // Test edge cases
    {
        const Point start(-0.5, 1, 0);
        const auto params = test.world_to_params(start);
        const auto start_point = test.face_->evaluate(params);

        // Zero length path
        try {
            [[maybe_unused]] auto path = test.face_->create_path(
                start_point, Vector(1, 0, 0), 0.0);
            assert(false && "Expected exception for zero length");
        } catch (const std::invalid_argument&) {}

        // Perpendicular direction
        try {
            [[maybe_unused]] auto path = test.face_->create_path(
                start_point, Vector(0, 1, 0), 1.0);
            assert(false && "Expected exception for perpendicular direction");
        } catch (const std::runtime_error&) {}
    }
}

// Test path transitions between cube faces
void test_cube_face_transitions() {
    CubeTest test;
    auto front = test.cube_.get_surface(0);
    assert(front);
    
    // Test front to right transition
    {
        const Point start(0.5, 1, 0);
        const Vector dir = Vector(1, -0.5, 0).normalized();
        const double length = std::sqrt(1.25);

        const auto params = CubeFaceTest::world_to_params(start);
        const auto start_point = front->evaluate(params);
        auto path = test.cube_.create_path(start_point, dir, length);

        // Check transition point
        const auto trans_pt = path->evaluate(0.5);
        const Point expected_trans(1, 0.75, 0);
        assert(approx_equal(trans_pt.world_pos(), expected_trans));
    }

    // Test complete circuit
    {
        const Point start(0, 1, 0);
        const Vector dir(1, 0, 0);
        const double length = 8.0;  // Full circuit

        const auto params = CubeFaceTest::world_to_params(start);
        const auto start_point = front->evaluate(params);
        auto path = test.cube_.create_path(start_point, dir, length);

        // Check face transitions
        const std::vector<std::pair<double, Point>> transitions = {
            {1.0, Point(1, 0, 0)},    // Front to right
            {3.0, Point(0, -1, 0)},   // Right to back
            {5.0, Point(-1, 0, 0)},   // Back to left
            {7.0, Point(0, 1, 0)}     // Left to front
        };

        for (const auto& [t, expected] : transitions) {
            const auto pt = path->evaluate(t);
            assert(approx_equal(pt.world_pos(), expected));
        }
    }
}

} // namespace shap::test

int main() {
    std::cout << "Running path tests...\n\n";
    
    try {
        shap::test::test_cube_face_paths();
        std::cout << "\n";
        
        shap::test::test_cube_face_transitions();
        std::cout << "\n";
        
        std::cout << "All tests passed!\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << "\n";
        return 1;
    }
}
