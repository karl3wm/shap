# Path Length Invariants

---
status: OPEN
test_file: tests/path_length_tests.cpp
test_name: "test_path_length_invariants() - Test 1: Parameter Space Distance"
first_seen: 2025-01-22
last_updated: 2025-01-22
resolution: null
related_issues:
  - length_scaling.md
  - metric_tensor_analysis.md
  - space_transformations.md
---

## Core Problem

The test reveals that our path length is not being preserved. For a path of requested length L=1.0:
```
t = 0.2 should give distance = 0.2
but gives distance = 0.3
```

This indicates our parameter space transformation is not maintaining the correct relationship between path parameter t and world space distance.

## Mathematical Relationships

### 1. Surface Parameterization (P)
- Maps [0,1]×[0,1] → ℝ³
- For SquareFace: P(u,v) = origin + u*du + v*dv
- Where |du| = 2, |dv| = 2 in our test case

### 2. Path Parameter (t)
- Maps [0,1] → World Space Distance
- Must maintain: |P(t) - P(0)| = t*L
- Where L is requested path length

### 3. Parameter Space Velocity
For a path moving in the u direction:
```
du/dt = 1/(2|du|)  // Factor of 2 because du has length 2
```

This means:
- To move 1 unit in world space
- We need to move 0.5 units in parameter space
- Because du has length 2

## Current Implementation Issue

We're currently computing:
```cpp
param_length = length / du_scale;  // Scale by du since moving in x direction
end_params = start_params + du_param * param_length;
```

But this doesn't account for the fact that:
1. Parameter derivatives (du_param, dv_param) are already scaled by the metric
2. The parameter space needs to maintain constant world space velocity

## Required Invariants

1. Parameter Space Distance
```
Δu = L/(2|du|)  // For motion in u direction
```

2. World Space Distance
```
|P(t) - P(0)| = t*L
```

3. Constant Speed
```
d/dt|P(t) - P(0)| = L
```

## Solution Requirements

1. Parameter Derivatives
- Must account for surface metric
- Must preserve requested length
- Must maintain constant speed

2. Path Evaluation
- Must verify distance at each t
- Must reach exactly requested length
- Must maintain constant speed

3. Testing
- Must verify these invariants for:
  * Different scale factors
  * Different directions
  * Different path lengths
