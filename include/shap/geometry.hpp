#pragma once
#include <cmath>
#include <array>

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

// 2x2 matrix for metric tensor
struct MetricTensor {
    double g11, g12, g21, g22;
    
    constexpr MetricTensor(double g11, double g12, double g21, double g22)
        : g11(g11), g12(g12), g21(g21), g22(g22) {}
        
    constexpr double determinant() const {
        return g11 * g22 - g12 * g21;
    }
};

} // namespace shap
