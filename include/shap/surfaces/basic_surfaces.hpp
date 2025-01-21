#pragma once
#include "../surface.hpp"
#include "../surface_collection.hpp"
#include <cmath>

namespace shap::surfaces {

// Create a sphere surface using the function-based interface
inline std::shared_ptr<Surface> create_sphere(double radius = 1.0) {
    return Surface::create_with_all_derivatives(
        // Position function
        [radius](double u, double v) {
            return Point(
                radius * std::cos(u) * std::sin(v),
                radius * std::sin(u) * std::sin(v),
                radius * std::cos(v)
            );
        },
        // First derivative du
        [radius](double u, double v) {
            return Point(
                -radius * std::sin(u) * std::sin(v),
                radius * std::cos(u) * std::sin(v),
                0
            );
        },
        // First derivative dv
        [radius](double u, double v) {
            return Point(
                radius * std::cos(u) * std::cos(v),
                radius * std::sin(u) * std::cos(v),
                -radius * std::sin(v)
            );
        },
        // Second derivative duu
        [radius](double u, double v) {
            return Point(
                -radius * std::cos(u) * std::sin(v),
                -radius * std::sin(u) * std::sin(v),
                0
            );
        },
        // Second derivative duv
        [radius](double u, double v) {
            return Point(
                -radius * std::sin(u) * std::cos(v),
                radius * std::cos(u) * std::cos(v),
                0
            );
        },
        // Second derivative dvv
        [radius](double u, double v) {
            return Point(
                -radius * std::cos(u) * std::sin(v),
                -radius * std::sin(u) * std::sin(v),
                -radius * std::cos(v)
            );
        },
        Surface::SurfaceType::Smooth
    );
}

// Create a square face using the function-based interface
inline std::shared_ptr<Surface> create_square_face(
    std::function<Point(double,double)> transform
) {
    return Surface::create_with_derivatives(
        // Position function
        std::move(transform),
        // Numerical du derivative
        [transform](double u, double v) {
            const double h = 1e-7;
            auto p1 = transform(u + h, v);
            auto p2 = transform(u - h, v);
            return (p1 - p2) * (0.5 / h);
        },
        // Numerical dv derivative
        [transform](double u, double v) {
            const double h = 1e-7;
            auto p1 = transform(u, v + h);
            auto p2 = transform(u, v - h);
            return (p1 - p2) * (0.5 / h);
        },
        Surface::SurfaceType::Developable
    );
}

// Create a cube using the improved collection interface
inline SurfaceCollection create_cube(double size = 1.0) {
    SurfaceCollection cube;
    
    // Create faces with appropriate transformations
    cube.add(create_square_face(
        [size](double u, double v) {
            return Point(size * (2*u - 1), size, size * (2*v - 1));
        }), "front"
    ).add(create_square_face(
        [size](double u, double v) {
            return Point(size, size * (1 - 2*u), size * (2*v - 1));
        }), "right"
    ).add(create_square_face(
        [size](double u, double v) {
            return Point(size * (1 - 2*u), -size, size * (2*v - 1));
        }), "back"
    ).add(create_square_face(
        [size](double u, double v) {
            return Point(-size, size * (2*u - 1), size * (2*v - 1));
        }), "left"
    );
    
    // Connect faces with automatic edge detection
    cube.connect("front", "right")
        .along(Edge::Right)
        .with_type(ConnectionType::Linear)
        .build();
        
    cube.connect("right", "back")
        .along(Edge::Right)
        .with_type(ConnectionType::Linear)
        .build();
        
    cube.connect("back", "left")
        .along(Edge::Right)
        .with_type(ConnectionType::Linear)
        .build();
        
    cube.connect("left", "front")
        .along(Edge::Right)
        .with_type(ConnectionType::Linear)
        .build();
    
    return cube;
}

} // namespace shap::surfaces