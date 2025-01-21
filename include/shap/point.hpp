#pragma once
#include <cmath>

namespace shap {

// Basic 3D point/vector class
class Point {
public:
    double x, y, z;
    
    Point(double x = 0, double y = 0, double z = 0)
        : x(x), y(y), z(z) {}
    
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
    
    double length() const {
        return std::sqrt(dot(*this));
    }
    
    Point normalize() const {
        double len = length();
        if (len < 1e-10) return *this;
        return *this * (1.0 / len);
    }
};

// Alias for vectors (same as points)
using Vector = Point;

// 2D point for parameter space coordinates
struct Point2D {
    double x, y;
    Point2D(double x = 0, double y = 0) : x(x), y(y) {}
};

} // namespace shap