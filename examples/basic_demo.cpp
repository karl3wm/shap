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
                  << normal.y << ", " << normal.z << ")\n\n";
    }
}

int main() {
    std::cout << "=== Surface Creation and Path Generation Examples ===\n\n";
    
    // 1. Creating a simple sphere using function-based interface
    std::cout << "1. Function-Based Surface Creation:\n\n";
    
    auto sphere = surfaces::create_sphere(1.0);
    sphere->name = "sphere";
    
    // Evaluate some points on the sphere
    auto sphere_point = sphere->evaluate(0.0, M_PI/4);
    print_point_info(sphere_point, "Sphere point at (u=0, v=π/4)");
    
    // 2. Creating a custom surface with position function
    std::cout << "2. Custom Surface Creation:\n\n";
    
    auto torus = Surface::create(
        [](double u, double v) {
            const double R = 2.0; // major radius
            const double r = 0.5; // minor radius
            return Point(
                (R + r * std::cos(v)) * std::cos(u),
                (R + r * std::cos(v)) * std::sin(u),
                r * std::sin(v)
            );
        },
        Surface::SurfaceType::Smooth
    );
    torus->name = "torus";
    
    auto torus_point = torus->evaluate(0.0, 0.0);
    print_point_info(torus_point, "Torus point at (u=0, v=0)");
    
    // 3. Creating and connecting surfaces in a collection
    std::cout << "3. Surface Collection and Connections:\n\n";
    
    // Create a cube with automatic connections
    auto cube = surfaces::create_cube(1.0);
    
    // Create a path that transitions between faces
    auto start_point = SurfacePoint(
        "front",           // Surface name
        0.8, 0.5,         // Near right edge, middle height
        Point(0.8, 1.0, 0.0),  // Position on front face
        Vector(0, 1, 0),       // Normal (front face)
        Vector(1, 0, 0),       // du (along u)
        Vector(0, 0, 1)        // dv (along v)
    );
    
    // Create path that wraps around the cube
    auto path = cube.create_path(
        start_point,
        Vector(1, 0, 0),  // Move towards right edge
        4.0               // Long enough to cross multiple faces
    );
    
    std::cout << "Path transitioning across cube faces:\n";
    print_path_info(*path);
    
    return 0;
}
