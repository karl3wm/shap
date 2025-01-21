#pragma once
#include <cmath>

namespace shap {

// 3D point/vector type with geometric operations
class Point {
public:
    double x, y, z;
    
    Point(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
    
    Point operator+(const Point& other) const {
        return Point(x + other.x, y + other.y, z + other.z);
    }
    
    Point operator-(const Point& other) const {
        return Point(x - other.x, y - other.y, z - other.z);
    }
    
    Point operator*(double s) const {
        return Point(x * s, y * s, z * s);
    }
    
    double dot(const Point& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    
    Point cross(const Point& other) const {
        return Point(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }
    
    double norm_squared() const {
        return dot(*this);
    }
    
    double norm() const {
        return std::sqrt(norm_squared());
    }
    
    Point normalize() const {
        double n = norm();
        return n > 0 ? *this * (1.0 / n) : Point();
    }
};

using Vector = Point;

} // namespace shap