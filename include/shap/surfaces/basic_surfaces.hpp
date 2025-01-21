#pragma once
#include "../surface.hpp"
#include "../surface_collection.hpp"
#include <cmath>

namespace shap::surfaces {

// Helper to create a square face with transformation
class SquareFace : public Surface {
public:
    using TransformFunc = std::function<Point(double,double)>;
    
    SquareFace(TransformFunc transform) : transform_(std::move(transform)) {}
    
    Point evaluate_position(double u, double v) const override {
        return transform_(u, v);
    }
    
    Point evaluate_du(double u, double v) const override {
        const double h = 1e-7;
        return (evaluate_position(u + h, v) - evaluate_position(u - h, v)) * (0.5 / h);
    }
    
    Point evaluate_dv(double u, double v) const override {
        const double h = 1e-7;
        return (evaluate_position(u, v + h) - evaluate_position(u, v - h)) * (0.5 / h);
    }
    
    // Square faces are developable (can be flattened)
    SurfaceType surface_type() const override {
        return SurfaceType::Developable;
    }
    
    // UV coordinates are already in flattened space
    std::pair<double,double> develop_point(double u, double v) const override {
        return {u, v};
    }

private:
    TransformFunc transform_;
};

// Sphere surface with given radius
class Sphere : public Surface {
public:
    explicit Sphere(double radius = 1.0) : radius_(radius) {}
    
    Point evaluate_position(double u, double v) const override {
        // u: longitude [0, 2π]
        // v: latitude [0, π]
        return Point(
            radius_ * std::cos(u) * std::sin(v),
            radius_ * std::sin(u) * std::sin(v),
            radius_ * std::cos(v)
        );
    }
    
    Point evaluate_du(double u, double v) const override {
        return Point(
            -radius_ * std::sin(u) * std::sin(v),
            radius_ * std::cos(u) * std::sin(v),
            0
        );
    }
    
    Point evaluate_dv(double u, double v) const override {
        return Point(
            radius_ * std::cos(u) * std::cos(v),
            radius_ * std::sin(u) * std::cos(v),
            -radius_ * std::sin(v)
        );
    }

    // Second derivatives for exact metric computation
    Point evaluate_duu(double u, double v) const override {
        return Point(
            -radius_ * std::cos(u) * std::sin(v),
            -radius_ * std::sin(u) * std::sin(v),
            0
        );
    }

    Point evaluate_duv(double u, double v) const override {
        return Point(
            -radius_ * std::sin(u) * std::cos(v),
            radius_ * std::cos(u) * std::cos(v),
            0
        );
    }

    Point evaluate_dvv(double u, double v) const override {
        return Point(
            -radius_ * std::cos(u) * std::sin(v),
            -radius_ * std::sin(u) * std::sin(v),
            -radius_ * std::cos(v)
        );
    }

    bool has_second_derivatives() const override {
        return true;
    }
    
    // Sphere is a smooth surface
    SurfaceType surface_type() const override {
        return SurfaceType::Smooth;
    }

private:
    double radius_;
};

// Create a cube as a collection of connected faces
inline SurfaceCollection create_cube(double size = 1.0) {
    SurfaceCollection cube;
    
    // Create faces with appropriate transformations
    auto front = std::make_shared<SquareFace>(
        [size](double u, double v) {
            return Point(size * (2*u - 1), size, size * (2*v - 1));
        }
    );
    front->name = "front";
    
    auto right = std::make_shared<SquareFace>(
        [size](double u, double v) {
            return Point(size, size * (1 - 2*u), size * (2*v - 1));
        }
    );
    right->name = "right";
    
    auto back = std::make_shared<SquareFace>(
        [size](double u, double v) {
            return Point(size * (1 - 2*u), -size, size * (2*v - 1));
        }
    );
    back->name = "back";
    
    auto left = std::make_shared<SquareFace>(
        [size](double u, double v) {
            return Point(-size, size * (2*u - 1), size * (2*v - 1));
        }
    );
    left->name = "left";
    
    // Add faces to collection
    cube.add_surface(front, "front");
    cube.add_surface(right, "right");
    cube.add_surface(back, "back");
    cube.add_surface(left, "left");
    
    // Add face connections with proper handling of non-smooth transitions
    // Front -> Right
    cube.add_connection(
        front, right,
        [](const SurfacePoint& pt, const Vector& dir) {
            return pt.u >= 0.95 && dir.x > 0;
        },
        [](const SurfacePoint& pt) {
            return SurfacePoint(
                "right",            // Next surface
                0.0, pt.v,         // Map to left edge
                pt.position,       // Keep position
                Vector(1, 0, 0),   // Normal points right
                Vector(0, -1, 0),  // du points back
                Vector(0, 0, 1)    // dv points up
            );
        }
    );
    
    // Right -> Back
    cube.add_connection(
        right, back,
        [](const SurfacePoint& pt, const Vector& dir) {
            return pt.u >= 0.95 && dir.x < 0;
        },
        [](const SurfacePoint& pt) {
            return SurfacePoint(
                "back",            // Next surface
                0.0, pt.v,        // Map to left edge
                pt.position,      // Keep position
                Vector(0, -1, 0), // Normal points back
                Vector(-1, 0, 0), // du points left
                Vector(0, 0, 1)   // dv points up
            );
        }
    );
    
    // Back -> Left
    cube.add_connection(
        back, left,
        [](const SurfacePoint& pt, const Vector& dir) {
            return pt.u >= 0.95 && dir.x < 0;
        },
        [](const SurfacePoint& pt) {
            return SurfacePoint(
                "left",           // Next surface
                0.0, pt.v,       // Map to left edge
                pt.position,     // Keep position
                Vector(-1, 0, 0), // Normal points left
                Vector(0, 1, 0),  // du points front
                Vector(0, 0, 1)   // dv points up
            );
        }
    );
    
    // Left -> Front
    cube.add_connection(
        left, front,
        [](const SurfacePoint& pt, const Vector& dir) {
            return pt.u >= 0.95 && dir.x > 0;
        },
        [](const SurfacePoint& pt) {
            return SurfacePoint(
                "front",          // Next surface
                0.0, pt.v,       // Map to left edge
                pt.position,     // Keep position
                Vector(0, 1, 0),  // Normal points front
                Vector(1, 0, 0),  // du points right
                Vector(0, 0, 1)   // dv points up
            );
        }
    );
    
    return cube;
}

} // namespace shap::surfaces