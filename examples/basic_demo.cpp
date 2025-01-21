#include <shap/surface.hpp>
#include <iostream>
#include <iomanip>

using namespace shap;

void print_point_info(const SurfacePoint& pt, const char* label = nullptr) {
    if (label) std::cout << label << ":\n";
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "Parameters: (u=" << pt.u << ", v=" << pt.v << ")\n";
    std::cout << "Position: (" << pt.position.x << ", " 
              << pt.position.y << ", " << pt.position.z << ")\n";
    std::cout << "Normal: (" << pt.normal.x << ", " 
              << pt.normal.y << ", " << pt.normal.z << ")\n\n";
}

void print_path_info(const SurfacePath& path, const char* label = nullptr) {
    if (label) std::cout << label << ":\n";
    
    // Sample points along path
    for (double t = 0; t <= 1.0; t += 0.1) {
        std::cout << "t = " << t << ":\n";
        auto pt = path.evaluate(t);
        std::cout << "  Position: (" << pt.position.x << ", "
                  << pt.position.y << ", " << pt.position.z << ")\n";
        
        auto tan = path.tangent(t);
        std::cout << "  Tangent: (" << tan.x << ", "
                  << tan.y << ", " << tan.z << ")\n\n";
    }
}

int main() {
    // Create a cube with connected faces
    auto cube = surfaces::cube(1.0);
    
    // Create a path that follows the cube surface
    // Start on front face, near top right corner
    auto start_point = SurfacePoint(
        0.8, 0.8,  // Near top right corner in uv space
        Point(0.8, 1.0, 0.8),  // Position on front face
        Vector(0, 1, 0),       // Normal (front face)
        Vector(1, 0, 0),       // du (along u)
        Vector(0, 0, 1)        // dv (along v)
    );
    
    // Create path in positive x direction (will transition to right face)
    auto path = cube.create_path(
        start_point,
        Vector(1, 0, 0),  // Move in +x direction
        4.0               // Long enough to wrap around
    );
    
    std::cout << "=== Original Path ===\n\n";
    print_path_info(*path, "Base path");
    
    // Create offset path raised above surface
    auto offset_path = path->offset(0.2);
    
    std::cout << "=== Offset Path ===\n\n";
    print_path_info(*offset_path, "Raised path");
    
    // Create smoothed version with rounded corners
    auto smooth_path = offset_path->smooth(0.1);
    
    std::cout << "=== Smoothed Path ===\n\n";
    print_path_info(*smooth_path, "Final ribbon path");
    
    return 0;
}
