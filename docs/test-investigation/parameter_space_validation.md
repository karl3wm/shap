# Parameter Space Validation and Contracts

## Surface Parameter Space

### Contracts

1. Parameter Domain
   - All surface parameters (u,v) must be in [0,1]×[0,1]
   - This is enforced by Surface::validate_parameters()
   - Represents a normalized coordinate system regardless of surface scale

2. World Space Mapping
   - Each surface defines its own mapping from parameter space to world space
   - Must be continuous and differentiable within the parameter domain
   - Example: SquareFace maps [0,1]×[0,1] to a rectangle in 3D space
   ```cpp
   // SquareFace mapping:
   P(u,v) = origin + u*du + v*dv
   // where du,dv define the face size and orientation
   ```

3. Parameter Space Derivatives
   - du = ∂P/∂u: Rate of change in u direction
   - dv = ∂P/∂v: Rate of change in v direction
   - Must be non-zero and non-parallel (enforced in SquareFace constructor)

## Path Parameter Space

### Contracts

1. Path Parameter Domain
   - All path parameters t must be in [0,1]
   - t=0 corresponds to path start point
   - t=1 corresponds to path end point
   - Enforced by SurfacePath::validate_parameter()

2. Path Length Mapping
   - Given a requested world space length L
   - Path parameter t must map to world space distance d = t*L
   - This ensures uniform speed along the path

3. Path Continuity
   - Position must vary continuously with t
   - Tangent vector must vary continuously with t
   - Essential for smooth paths and transitions

### Validation Points

1. Parameter Space Entry
   - Surface::create_path must validate input parameters
   - Must check start point is on surface
   - Must check direction is not zero
   - Must check length is positive

2. Parameter Space Exit
   - Path::evaluate must validate t in [0,1]
   - Must return point on surface
   - Must maintain requested world space length

3. Coordinate Transformations
   - world_to_parameters must handle surface scale factors
   - Parameter derivatives must account for surface metric
   - Path sampling must preserve arc length

## Diagnostic Requirements

1. Path Creation
   ```cpp
   // Log key information during path creation:
   - Input parameters (start, direction, length)
   - Surface properties at start point
   - Parameter space direction computation
   - End point computation
   - Parameter space sampling
   ```

2. Path Evaluation
   ```cpp
   // Log for each evaluation:
   - Input parameter t
   - Parameter space coordinates (u,v)
   - World space position
   - Distance from start point
   - Comparison with expected position
   ```

3. Surface Properties
   ```cpp
   // Log at key points:
   - Surface scale factors (du, dv lengths)
   - Parameter space bounds
   - World space bounds
   - Metric tensor components
   ```

## Testing Requirements

1. Path Length Validation
   ```cpp
   // For any path with length L:
   - Distance(P(t), P(0)) ≈ t*L for all t
   - Distance(P(1), P(0)) ≈ L
   ```

2. Parameter Space Coverage
   ```cpp
   // Test paths that:
   - Stay within parameter bounds
   - Reach parameter bounds
   - Follow parameter lines (constant u or v)
   - Cross parameter lines (diagonal paths)
   ```

3. Scale Factor Tests
   ```cpp
   // Verify behavior with:
   - Unit scale factors
   - Non-unit scale factors
   - Unequal scale factors
   - Varying scale factors
