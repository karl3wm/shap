# Length Scaling in Parameter Space

---
status: OPEN
test_file: tests/path_tests.cpp, tests/path_length_tests.cpp, tests/space_transformation_tests.cpp
test_name: "Multiple related failures in path length scaling"
first_seen: 2025-01-22
last_updated: 2025-01-22
resolution: null
related_issues:
  - path_length_invariants.md
  - metric_tensor_analysis.md
  - space_transformations.md
---

## Problem Analysis

A systematic investigation across multiple test files reveals a fundamental issue with length scaling between world space and parameter space. The issue manifests in path creation/evaluation, while basic coordinate transformations work correctly.

### Test Progression

1. path_tests.cpp: Initial failure showing paths end at wrong positions
2. path_length_tests.cpp: Confirms parameter space distances are wrong (0.75 vs 0.25)
3. space_transformation_tests.cpp: Verifies issue is specific to path creation/evaluation
   - Tests 1&2 PASS: Basic parameter<->world mappings work correctly
   - Test 3 FAILS: Path distances are 1.5x too large at t=0.2


The test failure reveals a fundamental issue with length scaling between world space and parameter space. Here's a concrete example from the failing test:

### Surface Setup
```
SquareFace:
- Parameter space: [0,1]×[0,1]
- World space: 2×2 unit square
- du = (2,0,0): u=1 maps to 2 units in x
- dv = (0,0,2): v=1 maps to 2 units in z
```

### Test Case
```
Start point:
- World: (-0.5, 1, 0)
- Parameters: (u=0.25, v=0.5)

Desired movement:
- Direction: (1,0,0) [x direction]
- Length: 1.0 units

Expected result:
- End at (0.5, 1, 0)
- Change in x = 1.0
- Required change in u = 0.5 [since du length = 2.0]
- Target u = 0.75

Actual result:
- Ends at (1.0, 1, 0)
- Change in x = 1.5
- Actual change in u = 0.75
```

### Key Insight
When creating paths on scaled surfaces:
1. World space lengths must be converted to parameter space
2. For a surface vector of length L:
   - Moving 1 unit in world space
   - Requires parameter change of 1/L
3. Must account for surface scale in both:
   - Initial direction computation
   - Path length computation

This explains why our current implementation overshoots - it's not properly accounting for the surface scale factors when converting between world space and parameter space lengths.
