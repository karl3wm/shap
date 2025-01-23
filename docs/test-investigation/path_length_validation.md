# Path Length Validation

This document consolidates key insights about path length calculations and validation from test investigations.

## Current Test Status

1. Passing Tests
- Parameter to World Space Mapping
  * Maps (0,0) → (-1,1,-1)
  * Maps (1,0) → (1,1,-1)
  * Maps (0,1) → (-1,1,1)
  * Maps (0.5,0.5) → (0,1,0)
- World to Parameter Space Mapping
  * All inverse mappings work correctly
  * Confirms basic coordinate transformations are sound

2. Failing Tests
- Path Length Preservation
  * Start: (-0.5,1,0)
  * Direction: (1,0,0)
  * Length: 1.0
  * Issue: Distance not preserved during path evaluation
  * Metric tensor values show correct computation:
    - g_uu = 4.0
    - g_uv = 0.0
    - g_vv = 4.0
    - det(g) = 16.0

## Core Requirements

1. Parameter Space to World Space Consistency
- Basic coordinate mappings are verified working
- Surface parameterization is mathematically correct
- Metric tensor computation is correct (g_uu=4, g_uv=0, g_vv=4)

2. Path Length Preservation (Currently Failing)
- For a path of requested length L=1.0:
  * |P(t) - P(0)| = t*L must hold for all t ∈ [0,1]
  * Current failure: at t=0.2, distance=0.3 (should be 0.2)
  * Parameter derivatives show correct metric scaling:
    - du/dt = 0.25 (accounts for |du|=2)
    - dv/dt = 0.0 (correct for x-direction motion)
- Scale factors appear to be applied correctly in metric tensor
- Issue may be in path evaluation rather than metric computation

## Test Cases (From Current Output)

1. Basic Coordinate Transformations (All Passing)
```cpp
// Parameter to World (Verified)
(0.0, 0.0) -> (-1.0, 1.0, -1.0)
(1.0, 0.0) -> (1.0, 1.0, -1.0)
(0.0, 1.0) -> (-1.0, 1.0, 1.0)
(0.5, 0.5) -> (0.0, 1.0, 0.0)

// World to Parameter (Verified)
(-1.0, 1.0, -1.0) -> (0.0, 0.0)
(1.0, 1.0, -1.0) -> (1.0, 0.0)
(-1.0, 1.0, 1.0) -> (0.0, 1.0)
(0.0, 1.0, 0.0) -> (0.5, 0.5)
```

2. Path Length Verification (Currently Failing)
```cpp
Path Configuration:
- Start: (-0.5, 1.0, 0.0)
- Direction: (1.0, 0.0, 0.0)
- Length: 1.0

Current Behavior (from test output):
- At t = 0.2: distance = 0.3 (should be 0.2)
- Parameter space coordinates correct:
  * Start: u=0.25, v=0.5
  * End: u=0.5, v=0.5
  * Delta: du=0.25, dv=0.0
- Path sampling shows linear progression
- Distance accumulation appears to be the issue
```

3. Scale Factor Handling
```cpp
For surface with |u_basis| = 2:
- Moving 1 unit in world space
- Requires Δu = 1/2 in parameter space
- Must account for metric tensor
```

## Implementation Requirements

1. Path Creation
- Use metric tensor for parameter derivatives
- Properly normalize directions
- Account for surface scale without double-counting

2. Path Evaluation
- Maintain constant world space velocity
- Preserve requested path length
- Handle scale factors consistently

3. Testing Strategy
- Verify coordinate transformations
- Check path lengths at multiple t values
- Test with different:
  * Scale factors
  * Directions
  * Path lengths

## Related Test Files
- tests/path_tests.cpp
- tests/path_length_tests.cpp
- tests/space_transformation_tests.cpp
