#include <shap/surface.hpp>
#include <iostream>
#include <iomanip>

using namespace shap;

void print_surface_info(const auto& surface, double u, double v) {
    std::cout << std::fixed << std::setprecision(6);
    
    auto p = surface(u, v);
    std::cout << "Point: (" << p.x << ", " << p.y << ", " << p.z << ")\n";
    
    auto n = surface.normal(u, v);
    std::cout << "Normal: (" << n.x << ", " << n.y << ", " << n.z << ")\n";
    
    auto metric = surface.metric_tensor(u, v);
    std::cout << "Metric tensor:\n";
    std::cout << "  [" << metric.g11 << " " << metric.g12 << "]\n";
    std::cout << "  [" << metric.g21 << " " << metric.g22 << "]\n";
    
    auto K = surface.gaussian_curvature(u, v);
    std::cout << "Gaussian curvature: " << K << "\n\n";
}

int main() {
    // Create a unit sphere
    constexpr auto sphere = surfaces::sphere(1.0);
    
    std::cout << "=== Sphere Analysis ===\n";
    // Check equator point
    std::cout << "At (u,v) = (0.000, 1.571) [equator point]:\n";
    print_surface_info(sphere, 0.0, M_PI/2);
    
    // Check north pole
    std::cout << "At (u,v) = (0.000, 0.000) [north pole]:\n";
    print_surface_info(sphere, 0.0, 0.0);
    
    // Create a cube
    constexpr auto cube = surfaces::cube(1.0);
    
    std::cout << "=== Cube Analysis ===\n";
    // Check center of front face
    std::cout << "At (u,v) = (0.000, 0.500) [front face center]:\n";
    print_surface_info(cube, 0.0, 0.5);
    
    // Check edge between front and right faces
    std::cout << "At (u,v) = (0.167, 0.500) [front-right edge]:\n";
    print_surface_info(cube, 1.0/6.0, 0.5);
    
    return 0;
}
