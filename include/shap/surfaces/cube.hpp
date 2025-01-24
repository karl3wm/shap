#pragma once
#include "shap/coord.hpp"
#include "shap/surfaces/flat_patch.hpp"
#include "shap/surface_collection.hpp"
#include <array>
#include <string_view>

namespace shap {
namespace surfaces {

namespace detail {
    // Face parameters for cube construction
    struct FaceParams {
        std::string_view name;  // Use string_view for compile-time strings
        Point origin;          // Origin at corner
        Vector du;             // Edge vector for u direction
        Vector dv;             // Edge vector for v direction
    };

    // Helper to create face parameters
    [[nodiscard]] constexpr FaceParams make_face_params(
        std::string_view name,
        const Point& origin,
        const Vector& du,
        const Vector& dv
    ) noexcept {
        return FaceParams{name, origin, du, dv};
    }

    // Helper to connect faces
    inline void connect_faces(
        SurfaceCollection& cube,
        Surface3D* face1,
        Surface3D* face2,
        ParamIndex param1,
        ParamBound bound1,
        ParamIndex param2,
        ParamBound bound2,
        int orientation
    ) {
        EdgeDescriptor edge1{param1, bound1, 0.0};
        EdgeDescriptor edge2{param2, bound2, 0.0};
        cube.add_connection(face1, edge1, face2, edge2, orientation);
    }
} // namespace detail

/**
 * Create a cube centered at the origin with given size.
 *
 * The cube is composed of six flat patches connected along their edges.
 * Each face is parameterized in its own [0,1]×[0,1] domain.
 *
 * Face Layout:
 * - Front:  y = +size
 * - Right:  x = +size
 * - Back:   y = -size
 * - Left:   x = -size
 * - Top:    z = -size
 * - Bottom: z = +size
 *
 * @param size Half-length of cube edges (must be positive)
 * @return Surface collection representing the cube
 * @throws std::invalid_argument if size <= 0
 */
[[nodiscard]] inline SurfaceCollection create_cube(double size = 1.0) {
    if (size <= 0) {
        throw std::invalid_argument("Cube size must be positive");
    }

    SurfaceCollection cube;
    const double double_size = 2.0 * size;
    
    // Define face parameters
    constexpr size_t NUM_FACES = 6;
    const std::array<detail::FaceParams, NUM_FACES> faces{{
        // Front face (y = size)
        detail::make_face_params(
            "front",
            Point(-size, size, -size),     // top-left corner
            Vector(double_size, 0, 0),      // u: left to right (+x)
            Vector(0, 0, double_size)       // v: top to bottom (+z)
        ),
        
        // Right face (x = size)
        detail::make_face_params(
            "right",
            Point(size, size, -size),      // top-left corner
            Vector(0, -double_size, 0),     // u: back (-y)
            Vector(0, 0, double_size)       // v: top to bottom (+z)
        ),
        
        // Back face (y = -size)
        detail::make_face_params(
            "back",
            Point(size, -size, -size),     // top-left corner
            Vector(-double_size, 0, 0),     // u: left (-x)
            Vector(0, 0, double_size)       // v: top to bottom (+z)
        ),
        
        // Left face (x = -size)
        detail::make_face_params(
            "left",
            Point(-size, -size, -size),    // top-left corner
            Vector(0, double_size, 0),      // u: forward (+y)
            Vector(0, 0, double_size)       // v: top to bottom (+z)
        ),
        
        // Top face (z = -size)
        detail::make_face_params(
            "top",
            Point(-size, -size, -size),    // back-left corner
            Vector(double_size, 0, 0),      // u: right (+x)
            Vector(0, double_size, 0)       // v: forward (+y)
        ),
        
        // Bottom face (z = size)
        detail::make_face_params(
            "bottom",
            Point(-size, -size, size),     // back-left corner
            Vector(double_size, 0, 0),      // u: right (+x)
            Vector(0, double_size, 0)       // v: forward (+y)
        )
    }};
    
    // Create faces and store pointers
    std::array<Surface3D*, NUM_FACES> face_ptrs;
    for (size_t i = 0; i < faces.size(); ++i) {
        auto surface = create_flat_patch(
            faces[i].origin,
            faces[i].du,
            faces[i].dv
        );
        face_ptrs[i] = surface.get();
        cube.add_surface(std::move(surface));
    }
    
    // Connect side faces (cyclic connections)
    const std::array<size_t, 4> side_faces{0, 1, 2, 3}; // front, right, back, left
    for (size_t i = 0; i < side_faces.size(); ++i) {
        const size_t next = (i + 1) % side_faces.size();
        detail::connect_faces(
            cube,
            face_ptrs[side_faces[i]], face_ptrs[side_faces[next]],
            ParamIndex::U, ParamBound::Upper,
            ParamIndex::U, ParamBound::Lower,
            1
        );
    }
    
    // Connect top face (index 4)
    const std::array<std::pair<ParamIndex, ParamBound>, 4> top_connections{{
        {ParamIndex::V, ParamBound::Upper},  // front
        {ParamIndex::U, ParamBound::Upper},  // right
        {ParamIndex::V, ParamBound::Lower},  // back
        {ParamIndex::U, ParamBound::Lower}   // left
    }};
    
    for (size_t i = 0; i < side_faces.size(); ++i) {
        detail::connect_faces(
            cube,
            face_ptrs[4], face_ptrs[side_faces[i]], // 4 is top face
            top_connections[i].first, top_connections[i].second,
            ParamIndex::U, ParamBound::Lower,
            i < 2 ? 1 : -1
        );
    }
    
    // Connect bottom face (index 5)
    for (size_t i = 0; i < side_faces.size(); ++i) {
        detail::connect_faces(
            cube,
            face_ptrs[5], face_ptrs[side_faces[i]], // 5 is bottom face
            top_connections[i].first, top_connections[i].second,
            ParamIndex::V, ParamBound::Upper,
            i < 2 ? 1 : -1
        );
    }
    
    return cube;
}

} // namespace surfaces
} // namespace shap
