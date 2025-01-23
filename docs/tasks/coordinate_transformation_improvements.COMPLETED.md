# [COMPLETED 2025-01-23] Coordinate Transformation Improvements

Status: ✅ Completed
Date: January 23, 2025
Changes:
- Implemented robust world_to_local transformation
- Added comprehensive test suite for coordinate transformations
- Fixed handling of degenerate cases
- Improved error messages and validation
- All test cases passing

# Coordinate Transformation Improvements

## Overview
Improve the coordinate transformation interface to be more mathematically precise, better documented, and more consistent with the codebase's validation approach. Introduce clear naming conventions for different coordinate spaces.

## Current Issues
1. Implicit assumptions about points lying "exactly" on the surface
2. Separate handling of degenerate cases with custom epsilon values
3. Unclear documentation of Cramer's rule implementation
4. Inconsistent validation approach
5. Multiple functions with overlapping functionality
6. Ambiguous naming of coordinate spaces

## Migration Notes
- ParamPoint2 replaces ParameterCoordinate
- ParamPoint3 is a new type that combines parameter coordinates with normal distance
- GeometryPoint2 replaces both GeometricProperties and SurfacePoint
- WorldPoint3 and WorldVector3 replace Point and Vector respectively
- Each type moves to its own file following one-class-per-file norm
- Types will move out of types.hpp to their own files

## Files to Delete/Replace
After migration is complete, these files will be removed:
1. include/shap/point.hpp (replaced by WorldPoint3 and WorldVector3)
2. include/shap/surface_point.hpp (replaced by GeometryPoint2)
3. include/shap/types.hpp (contents moved to individual files)

## New Files
1. include/shap/world_point3.hpp
2. include/shap/world_vector3.hpp
3. include/shap/param_index.hpp
4. include/shap/param_bound.hpp
5. include/shap/edge_descriptor.hpp
6. include/shap/surface_type.hpp

## Epsilon Usage

Two distinct types of epsilon checks in coordinate transformations:

1. **Validation Epsilons** (Use ValidationConfig)
   - Checking if basis vectors are nearly parallel (throws exception)
   - Checking if a point lies exactly on surface
   These represent invariants that must hold for valid operation.

2. **Algorithm Epsilons** (Keep as parameters)
   - Parameter bound checks in path solving
   - Degenerate case handling where alternative computation paths exist
   These represent tunable tolerances in algorithms that don't throw exceptions.

Example in FlatPatch:
```cpp
// Validation epsilon - throws on invalid state
if (det < ValidationConfig::instance().vector_length_epsilon()) {
    throw std::invalid_argument("basis vectors are nearly parallel");
}

// Algorithm epsilon - handles degenerate case
if (std::abs(d_param) > parameter_bound_epsilon) {
    // Regular path
} else {
    // Alternative computation
}
```

## Coordinate Space Design

### 1. Point Types
Replace current coordinate representations with clear, strongly-typed points:

```cpp
// Points in local parameter space
class ParamPoint2 {
    double u, v;  // [0,1] × [0,1] coordinates
public:
    explicit ParamPoint2(double u, double v);  // With validation
    [[nodiscard]] double u() const noexcept { return u_; }
    [[nodiscard]] double v() const noexcept { return v_; }
};

// Points in local parameter space with normal offset
class ParamPoint3 {
    double u, v;      // [0,1] × [0,1] coordinates
    double normal;    // Signed normal distance
public:
    explicit ParamPoint3(double u, double v, double normal);
    [[nodiscard]] ParamPoint2 to_r2() const;  // Project to surface
    [[nodiscard]] double u() const noexcept { return u_; }
    [[nodiscard]] double v() const noexcept { return v_; }
    [[nodiscard]] double normal() const noexcept { return normal_; }
};

// Points with geometric properties
class GeometryPoint2 {
    ParamPoint2 local_pos;
    WorldPoint3 world_pos;       // Position in world space
    WorldVector3 world_normal;   // Surface normal vector
    WorldVector3 world_du;       // Derivative vector in u direction
    WorldVector3 world_dv;       // Derivative vector in v direction
public:
    GeometryPoint2(
        ParamPoint2 local,
        WorldPoint3 pos,
        WorldPoint3 normal,
        WorldPoint3 du,
        WorldPoint3 dv
    );
};
```

### 2. Surface Interface
Update the Surface class with clear coordinate transformation methods:

```cpp
class Surface {
public:
    /**
     * Convert a world space position to local coordinates.
     * 
     * This function computes three coordinates that fully describe a point's position
     * relative to the surface:
     * - u,v: Param parameter coordinates in [0,1]×[0,1]
     * - normal: Signed distance along surface normal vector
     *
     * For points on the surface, normal will be 0 (within ValidationConfig::vector_length_epsilon).
     * Positive normal indicates the point is on the positive side of the surface
     * (in the direction of the normal vector).
     *
     * @param pos World space position to convert
     * @return ParamPoint3 containing local coordinates
     * @throws std::invalid_argument if coordinate computation fails
     */
    [[nodiscard]] virtual ParamPoint3 world_to_local(const WorldPoint3& pos) const = 0;
    
    /**
     * Convert a world space position to surface parameter coordinates.
     * Projects the point onto the surface along the normal direction.
     *
     * @param pos World space position to convert
     * @return ParamPoint2 containing parameter coordinates
     * @throws std::invalid_argument if coordinate computation fails
     */
    [[nodiscard]] virtual ParamPoint2 world_to_local_r2(const WorldPoint3& pos) const {
        return world_to_local(pos).to_r2();
    }
};
```

### 3. FlatPatch Implementation
Update FlatPatch with the new interface:

```cpp
class FlatPatch {
    [[nodiscard]] ParamPoint3 world_to_local(const WorldPoint3& pos) const override {
        const Vector rel_pos = pos - origin_;
        
        // Project point onto surface normal to get signed distance
        const double normal_dist = dot(rel_pos, normal_);
        
        // Project point onto surface plane
        const Vector planar_pos = rel_pos - normal_dist * normal_;
        
        // Solve for u,v using Cramer's rule
        const double det = cross(world_u_, world_v_).length();
        if (det < ValidationConfig::instance().vector_length_epsilon()) {
            throw std::invalid_argument(
                "Cannot compute local coordinates: basis vectors are nearly parallel"
            );
        }
        
        return ParamPoint3(
            dot(cross(planar_pos, world_v_), normal_) / det,  // u coordinate
            dot(cross(world_u_, planar_pos), normal_) / det,  // v coordinate
            normal_dist  // signed distance from surface
        );
    }
};
```

## Implementation Order
1. First implement and test the point types (ParamPoint2, ParamPoint3)
2. Then implement GeometryPoint2 as it depends on the other point types
3. Update Surface interface
4. Finally update FlatPatch implementation
This order minimizes dependency issues during the transition.

## Testing Strategy
1. Start with unit tests for the new point types
2. Add conversion tests between old and new types during transition
3. Add the comprehensive coordinate transformation tests
4. Finally add integration tests with path solving

## Test Cases
Add comprehensive tests in space_transformation_tests.cpp:

```cpp
TEST_CASE("FlatPatch coordinate transformation") {
    // 1. Basic orthogonal basis
    SECTION("orthogonal basis") {
        // Points exactly on surface
        auto local = patch.world_to_local(on_surface_point);
        CHECK(local.normal() == Approx(0.0));
        CHECK(local.u() == Approx(0.5));
        CHECK(local.v() == Approx(0.5));
        
        // Points above and below surface
        auto above = patch.world_to_local(offset_point);
        CHECK(above.normal() > 0);
        
        // Points outside parameter range
        auto outside = patch.world_to_local(far_point);
        CHECK(outside.u() > 1.0);
    }
    
    // 2. Non-orthogonal basis
    SECTION("skewed basis") {
        // Verify Cramer's rule handles skewed coordinates
        auto local = skewed_patch.world_to_local(test_point);
        CHECK(local.u() == Approx(expected_u));
        CHECK(local.v() == Approx(expected_v));
    }
    
    // 3. Edge cases
    SECTION("edge cases") {
        // Points very far from surface
        // Points very close to surface
        // Points at parameter space boundaries
    }
    
    // 4. Error cases
    SECTION("degenerate cases") {
        // Nearly parallel basis vectors
        CHECK_THROWS_AS(
            bad_patch.world_to_local(test_point),
            std::invalid_argument
        );
    }
}
```

## Benefits
1. Clear separation of coordinate spaces through strong typing
2. Consistent naming convention (local_/world_ prefixes)
3. Explicit dimensionality (R2/R3 suffixes)
4. Better support for common use cases:
   - Ray intersection calculations
   - Geometry alignment
   - Surface offset operations
   - Distance measurements
   - Extrusion planning

## Future Extensions
1. Support for higher dimensional parameter spaces (R4, etc.)
2. Volume parameterizations (ParamPoint3 → WorldPoint3)
3. Additional geometric properties as needed
4. Utility functions for common operations (projection, offset, etc.)
