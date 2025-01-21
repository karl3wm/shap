#pragma once
#include <cmath>
#include <array>
#include <functional>
#include <memory>
#include <vector>

namespace shap {

// Basic 3D point/vector type with constexpr operations
struct Point {
    double x, y, z;
    
    constexpr Point(double x = 0, double y = 0, double z = 0) 
        : x(x), y(y), z(z) {}
    
    constexpr Point operator+(const Point& other) const {
        return Point(x + other.x, y + other.y, z + other.z);
    }
    
    constexpr Point operator-(const Point& other) const {
        return Point(x - other.x, y - other.y, z - other.z);
    }
    
    constexpr Point operator*(double s) const {
        return Point(x * s, y * s, z * s);
    }
    
    constexpr double dot(const Point& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    
    constexpr Point cross(const Point& other) const {
        return Point(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }
    
    constexpr double norm_squared() const {
        return dot(*this);
    }
    
    constexpr Point normalize() const {
        double n = std::sqrt(norm_squared());
        return *this * (1.0 / n);
    }
};

using Vector = Point;

// Surface point with parameter coordinates and geometric data
struct SurfacePoint {
    std::string surface_name;  // Name of containing surface
    double u, v;              // Parameter coordinates
    Point position;           // 3D position
    Vector normal;            // Surface normal
    Vector du, dv;           // Tangent vectors
    
    // Construct from parameters and geometric data
    SurfacePoint(
        std::string surface,
        double u, double v,
        Point pos,
        Vector n,
        Vector du, Vector dv
    ) : surface_name(std::move(surface)),
        u(u), v(v),
        position(pos),
        normal(n),
        du(du), dv(dv) {}
        
    // Get tangent vector in given parameter direction
    Vector tangent(double du_component, double dv_component) const {
        return (du * du_component + dv * dv_component).normalize();
    }
};

// 2x2 matrix for metric tensor with Riemannian geometry operations
struct MetricTensor {
    double g11, g12, g21, g22;
    
    constexpr MetricTensor(double g11, double g12, double g21, double g22)
        : g11(g11), g12(g12), g21(g21), g22(g22) {}
    
    // Compute Christoffel symbols for geodesic equations
    std::array<double,2> christoffel_first(int i, int j, int k) const;
    std::array<double,2> christoffel_second(int i) const;
    
    // Convert tangent vector components between coordinate systems
    std::pair<double,double> raise_indices(double v1, double v2) const {
        double det = determinant();
        return {
            (g22 * v1 - g12 * v2) / det,
            (-g21 * v1 + g11 * v2) / det
        };
    }
    
    constexpr double determinant() const {
        return g11 * g22 - g12 * g21;
    }
};

// Abstract base class for parametric paths on surfaces
class SurfacePath {
public:
    virtual ~SurfacePath() = default;
    
    // Evaluate path at parameter t in [0,1]
    virtual SurfacePoint evaluate(double t) const = 0;
    
    // Get tangent vector at parameter t
    virtual Vector tangent(double t) const = 0;
    
    // Create offset path at constant distance
    virtual std::unique_ptr<SurfacePath> offset(double distance) const = 0;
    
    // Create path with smoothed corners using circular arcs
    virtual std::unique_ptr<SurfacePath> smooth(double radius) const = 0;
};

// Geodesic path between two points on a surface
class GeodesicPath : public SurfacePath {
    // Implementation will use Riemannian connection to compute geodesics
};

// Path following constant direction on surface
class DirectionalPath : public SurfacePath {
    // Implementation will use parallel transport to maintain direction
};

// Composite path made up of multiple segments
class CompositePath : public SurfacePath {
    std::vector<std::unique_ptr<SurfacePath>> segments;
public:
    void add_segment(std::unique_ptr<SurfacePath> segment) {
        segments.push_back(std::move(segment));
    }
};

// Helper to create paths between surface points
std::unique_ptr<SurfacePath> create_geodesic_path(
    const SurfacePoint& start,
    const SurfacePoint& end
);

// Helper to create paths in constant direction
std::unique_ptr<SurfacePath> create_directional_path(
    const SurfacePoint& start,
    const Vector& direction,
    double length
);

} // namespace shap
