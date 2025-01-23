#include "shap/coord.hpp"
#include <shap/surfaces/basic_surfaces.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace shap::demo {

// Print utilities
namespace {
    constexpr int PRECISION = 6;
    constexpr double SAMPLE_INTERVAL = 0.1;
    
    // Helper to print section headers
    void print_section(std::string_view title) {
        std::cout << "\n=== " << title << " ===\n\n";
    }

    // Helper to print point information
    void print_point_info(const SurfacePoint& pt, std::string_view label = "") {
        if (!label.empty()) std::cout << label << ":\n";
        std::cout << std::fixed << std::setprecision(PRECISION);
        
        std::cout << "Surface: " << (pt.surface() ? "valid" : "null") << "\n"
                 << "Parameters: (u=" << pt.parameter_coordinates().u() 
                 << ", v=" << pt.parameter_coordinates().v() << ")\n"
                 << "Position: (" << pt.world_pos().x << ", " 
                 << pt.world_pos().y << ", " << pt.world_pos().z << ")\n"
                 << "Normal: (" << pt.world_normal().x << ", " 
                 << pt.world_normal().y << ", " << pt.world_normal().z << ")\n";
        
        if (auto edge = pt.get_edge_descriptor()) {
            std::cout << "Edge: param=" << static_cast<int>(edge->param)
                     << " bound=" << static_cast<int>(edge->bound)
                     << " value=" << edge->edge_param << "\n";
        }
        std::cout << "\n";
    }

    // Helper to print path information
    void print_path_info(const SurfacePath& path, std::string_view label = "") {
        if (!label.empty()) std::cout << label << ":\n";
        
        // Pre-compute sample points for efficiency
        std::vector<std::pair<double, SurfacePoint>> samples;
        samples.reserve(static_cast<size_t>(1.0 / SAMPLE_INTERVAL) + 1);
        
        for (double t = 0; t <= 1.0; t += SAMPLE_INTERVAL) {
            samples.emplace_back(t, path.evaluate(t));
        }
        
        // Print sample points
        for (const auto& [t, pt] : samples) {
            const auto normal = path.normal(t);
            
            std::cout << "t = " << std::fixed << std::setprecision(3) << t << ":\n"
                     << "  Surface: " << (pt.surface() ? "valid" : "null") << "\n"
                     << "  Position: " << std::fixed << std::setprecision(PRECISION)
                     << "(" << pt.world_pos().x << ", " << pt.world_pos().y << ", " 
                     << pt.world_pos().z << ")\n"
                     << "  Normal: (" << normal.x << ", " << normal.y << ", " 
                     << normal.z << ")\n";
            
            if (auto edge = pt.get_edge_descriptor()) {
                std::cout << "  Edge: param=" << static_cast<int>(edge->param)
                         << " bound=" << static_cast<int>(edge->bound)
                         << " value=" << edge->edge_param << "\n";
            }
            std::cout << "\n";
        }
    }

    // Helper to measure execution time
    template<typename F>
    double time_operation(F&& func) {
        const auto start = std::chrono::high_resolution_clock::now();
        func();
        const auto end = std::chrono::high_resolution_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        return static_cast<double>(duration.count()) / 1000.0;  // Convert to milliseconds
    }

    // Example 1: Basic sphere creation and evaluation
    void demo_basic_sphere() {
        print_section("Basic Sphere Creation and Evaluation");
        
        auto sphere = surfaces::create_sphere(1.0);
        
        // Evaluate points at different locations
        const std::vector<std::pair<double, double>> sample_points = {
            {0.0, M_PI/4},    // Top quadrant
            {M_PI/2, M_PI/2}, // Equator
            {M_PI, 3*M_PI/4}  // Bottom quadrant
        };
        
        for (const auto& [u, v] : sample_points) {
            auto point = sphere->evaluate(ParameterCoordinate{u/(2*M_PI), v/M_PI});
            print_point_info(point, "Sphere point at (u=" + 
                std::to_string(u) + ", v=" + std::to_string(v) + ")");
        }
    }

    // Example 2: Custom torus surface
    void demo_custom_torus() {
        print_section("Custom Torus Surface Creation");
        
        constexpr double R = 2.0; // major radius
        constexpr double r = 0.5; // minor radius
        
        // Position function
        auto position_func = [R, r](const ParameterCoordinate& param) {
            const double u = param.u() * 2 * M_PI;
            const double v = param.v() * 2 * M_PI;
            return Point(
                (R + r * std::cos(v)) * std::cos(u),
                (R + r * std::cos(v)) * std::sin(u),
                r * std::sin(v)
            );
        };
        
        // Analytical first derivatives
        auto du_func = [R, r](const ParameterCoordinate& param) {
            const double u = param.u() * 2 * M_PI;
            const double v = param.v() * 2 * M_PI;
            return Point(
                -(R + r * std::cos(v)) * std::sin(u),
                (R + r * std::cos(v)) * std::cos(u),
                0
            );
        };
        
        auto dv_func = [R, r](const ParameterCoordinate& param) {
            const double u = param.u() * 2 * M_PI;
            const double v = param.v() * 2 * M_PI;
            return Point(
                -r * std::sin(v) * std::cos(u),
                -r * std::sin(v) * std::sin(u),
                r * std::cos(v)
            );
        };
        
        // Analytical Gaussian curvature
        auto gaussian_curv = [R, r](const ParameterCoordinate& param) {
            const double v = param.v() * 2 * M_PI;
            const double cos_v = std::cos(v);
            return cos_v / ((R + r * cos_v) * r);
        };
        
        auto torus = Surface::create_with_derivatives(
            std::move(position_func),
            std::move(du_func),
            std::move(dv_func),
            std::move(gaussian_curv),
            std::nullopt,  // Mean curvature
            std::nullopt,  // Path solver
            SurfaceType::Smooth
        );
        
        // Sample points around the torus
        const std::vector<std::pair<double, double>> sample_points = {
            {0.0, 0.0},     // Outer equator
            {M_PI/2, 0.0},  // Side
            {0.0, M_PI/2}   // Top
        };
        
        for (const auto& [u, v] : sample_points) {
            auto point = torus->evaluate(ParameterCoordinate{u/(2*M_PI), v/(2*M_PI)});
            print_point_info(point, "Torus point at (u=" + 
                std::to_string(u) + ", v=" + std::to_string(v) + ")");
        }
    }

    // Example 3: Surface transitions on cube
    void demo_cube_transitions() {
        print_section("Surface Collection and Transitions");
        
        auto cube = surfaces::create_cube(1.0);
        
        auto front_face = cube.get_surface(0); // Front face is index 0
        if (!front_face) {
            throw std::runtime_error("Failed to get front face");
        }
        
        // Create paths with different starting points and directions
        const std::vector<std::tuple<double, double, Vector, double>> paths = {
            {0.9, 0.5, Vector(1, 0, 0), 2.0},    // Right transition
            {0.5, 0.1, Vector(0, 0, 1), 2.0},    // Bottom transition
            {0.1, 0.5, Vector(-1, 0, 0), 2.0}    // Left transition
        };
        
        for (const auto& [u, v, dir, length] : paths) {
            auto start_point = front_face->evaluate(ParameterCoordinate{u, v});
            print_point_info(start_point, "Starting point");
            
            auto path = cube.create_path(start_point, dir, length);
            print_path_info(*path, "Transition path");
        }
    }

    // Example 4: Performance measurements
    void demo_performance() {
        print_section("Performance Measurements");
        
        // Measure sphere creation and evaluation
        const double sphere_time = time_operation([]() {
            auto sphere = surfaces::create_sphere(1.0);
            for (int i = 0; i < 1000; ++i) {
                const double u = static_cast<double>(i) / 1000;
                const double v = static_cast<double>(i) / 1000;
                [[maybe_unused]] auto point = sphere->evaluate(ParameterCoordinate{u, v});
            }
        });
        std::cout << "Sphere creation and 1000 evaluations: " 
                 << sphere_time << "ms\n";
        
        // Measure cube path creation
        const double cube_time = time_operation([]() {
            auto cube = surfaces::create_cube(1.0);
            auto front = cube.get_surface(0); // Front face is index 0
            auto start = front->evaluate(ParameterCoordinate{0.5, 0.5});
            for (int i = 0; i < 100; ++i) {
                const double angle = 2 * M_PI * i / 100;
                const Vector dir(std::cos(angle), 0, std::sin(angle));
                [[maybe_unused]] auto path = cube.create_path(start, dir, 2.0);
            }
        });
        std::cout << "Cube creation and 100 paths: " << cube_time << "ms\n";
    }
}

} // namespace shap::demo

int main() {
    try {
        // Run all demos
        shap::demo::demo_basic_sphere();
        shap::demo::demo_custom_torus();
        shap::demo::demo_cube_transitions();
        shap::demo::demo_performance();
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
