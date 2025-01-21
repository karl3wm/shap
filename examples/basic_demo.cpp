#include <shap/surfaces/basic_surfaces.hpp>
#include <iostream>
#include <iomanip>

using namespace shap;

void print_point_info(const SurfacePoint& pt, const char* label = nullptr) {
    if (label) std::cout << label << ":\n";
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "Surface: " << pt.surface_name << "\n";
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
        auto normal = path.normal(t);
        std::cout << "  Surface: " << pt.surface_name << "\n";
        std::cout << "  Position: (" << pt.position.x << ", "
                  << pt.position.y << ", " << pt.position.z << ")\n";
        std::cout << "  Normal: (" << normal.x << ", "
                  << normal.y << ", " << normal.z << ")\n";
        auto tan = path.tangent(t);
        std::cout << "  Tangent: (" << tan.x << ", "
                  << tan.y << ", " << tan.z << ")\n\n";
    }
}

int main() {
    std::cout << "=== Surface Types and Geodesic Behavior ===\n\n";
    
    // Create surfaces with different geometric properties
    auto cube = surfaces::create_cube(1.0);
    auto sphere = std::make_shared<surfaces::Sphere>(1.0);
    sphere->name = "sphere";
    
    std::cout << "1. Cube Faces (Developable Surfaces):\n"
              << "   - Each face is flat and can be developed (flattened)\n"
              << "   - Geodesics within a face are straight lines\n"
              << "   - At edges, geodesics have direction changes\n"
              << "   - Transitions between faces are non-smooth\n\n";
    
    // Create path on cube with different transition types
    auto start_point = SurfacePoint(
        "front",           // Surface name
        0.8, 0.8,         // Near top right corner in uv space
        Point(0.8, 1.0, 0.8),  // Position on front face
        Vector(0, 1, 0),       // Normal (front face)
        Vector(1, 0, 0),       // du (along u)
        Vector(0, 0, 1)        // dv (along v)
    );
    
    std::cout << "2. Transition Types at Cube Edges:\n\n";
    
    try {
        std::cout << "Attempting geodesic transition (should fail):\n";
        auto path = cube.create_path(
            start_point,
            Vector(1, 0, 0),
            4.0
        );
        auto geodesic_path = path->offset(0.2, TransitionType::Geodesic);
        print_path_info(*geodesic_path, "Geodesic transitions (invalid for cube)");
    } catch (const std::runtime_error& e) {
        std::cout << "Error: " << e.what() << "\n"
                  << "This is expected because cube edges are non-smooth.\n\n";
    }
    
    std::cout << "Linear transition (simple corner):\n";
    auto path = cube.create_path(
        start_point,
        Vector(1, 0, 0),
        4.0
    );
    auto linear_path = path->offset(0.2, TransitionType::Linear);
    print_path_info(*linear_path, "Linear transitions");
    
    std::cout << "\nCircular arc transition (smooth corner):\n";
    auto circular_path = path->offset(0.2, TransitionType::Circular);
    print_path_info(*circular_path, "Circular arc transitions");
    
    // Demonstrate geodesics on smooth sphere
    std::cout << "\n3. Sphere (Smooth Surface):\n"
              << "   - Globally smooth surface\n"
              << "   - Geodesics follow great circles\n"
              << "   - No sharp transitions or edges\n\n";
    
    auto sphere_start = SurfacePoint(
        "sphere",
        0.0, M_PI/4,           // Parameters
        Point(0.7, 0, 0.7),    // Position
        Point(0.7, 0, 0.7),    // Normal
        Vector(0, 1, 0),       // du
        Vector(-0.7, 0, 0.7)   // dv
    );
    
    auto sphere_end = SurfacePoint(
        "sphere",
        M_PI/2, M_PI/4,       // Parameters
        Point(0, 0.7, 0.7),   // Position
        Point(0, 0.7, 0.7),   // Normal
        Vector(-1, 0, 0),     // du
        Vector(0, -0.7, 0.7)  // dv
    );
    
    try {
        auto sphere_geodesic = std::make_unique<GeodesicCurve>(
            sphere,
            sphere_start,
            sphere_end
        );
        print_path_info(*sphere_geodesic, "Sphere geodesic (great circle)");
    } catch (const std::runtime_error& e) {
        std::cout << "Error computing sphere geodesic: " << e.what() << "\n";
    }
    
    return 0;
}
