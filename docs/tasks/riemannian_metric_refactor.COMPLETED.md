# [COMPLETED 2025-01-23] Riemannian Metric Refactor

Status: ✅ Completed
Date: January 23, 2025
Changes:
- Moved function pointers into Surface as primary resources
- Created consolidated RiemannianMetric class
- Maintained mathematical documentation in code comments
- Preserved existing geometric operations
- Removed legacy metric.cpp/hpp files
- Updated build system to reflect changes
- Verified space transformation tests pass

## Objective
Refactor the metric-related functionality to:
1. Move function pointers into Surface as primary resources
2. Create a consolidated RiemannianMetric class
3. Maintain mathematical documentation in code comments
4. Preserve existing geometric operations

## Implementation Steps

1. Create include/shap/riemannian_metric.hpp:
   - Define RiemannianMetric class with constant components
   - Move relevant mathematical comments from metric.hpp
   - Include constructors for direct component specification and from Surface

2. Update include/shap/surface.hpp:
   - Add metric derivative function pointers
   - Document derivative calculations in constructor comments
   - Provide accessor methods for derivatives

3. Update existing surfaces:
   - Implement analytical derivatives for each surface type
   - Remove numerical approximations
   - Update constructors to provide derivative functions

4. Clean up metric.hpp/cpp:
   - Remove varying metric functionality
   - Move Christoffel symbol calculations to RiemannianMetric
   - Update documentation to reference new structure

## Code Structure

```cpp
// In riemannian_metric.hpp
class RiemannianMetric {
    // Components
    double g11, g12, g22;
    double dg11_du, dg11_dv;
    double dg12_du, dg12_dv;
    double dg22_du, dg22_dv;

    // Constructors
    RiemannianMetric(/* direct components */);
    explicit RiemannianMetric(const Surface& surface);

    // Core operations
    determinant();
    raise_indices();
    lower_indices();
    christoffel_first();
    christoffel_second();
    // etc.
};

// In surface.hpp
class Surface {
    // Function pointers
    DerivativeFunction du_fn_, dv_fn_;
    MetricDerivativeFunction dg11_du_fn_, dg11_dv_fn_;
    MetricDerivativeFunction dg12_du_fn_, dg12_dv_fn_;
    MetricDerivativeFunction dg22_du_fn_, dg22_dv_fn_;

    // Constructor and accessors
    Surface(/* function pointers */);
    du_at();
    dv_at();
    dg11_du();
    // etc.
};
```

## Testing
- Update existing tests to use new RiemannianMetric class
- Verify analytical derivatives match previous numerical results
- Ensure geodesic calculations remain accurate
- Test metric construction from different surface types

## Future Considerations
- Optimization opportunities by caching computed values
- Extension to other Riemannian geometry concepts
- Support for more complex surface types

## File Organization
Following the one-class-one-file norm:
1. include/shap/riemannian_metric.hpp - RiemannianMetric class
2. src/riemannian_metric.cpp - RiemannianMetric implementation
3. include/shap/surface.hpp - Updated Surface class
4. src/surface.cpp - Surface implementation
