# Path and Surface Unification

## Overview
Unifying paths and surfaces under a common manifold framework to leverage their shared characteristics:
- Parameter space to world space mapping
- Connection/transition handling
- Geometric information at points

## Current State
- Manifold base template implemented and working
- Surface implementations successfully converted to use Manifold framework
- Basic manifold operations verified through tests
- Coordinate system and geometric point templates in place
- Surface3D successfully inheriting from Manifold<2,3,WorldSpaceTag>

## Remaining Implementation

### 1. Path Implementation (In Progress)
Converting path classes to use Manifold framework:
```cpp
// Implemented structure
using WorldPath3D = Manifold<1, 3, WorldSpaceTag>;  // 1D -> 3D world
using ParamPath2D = Manifold<1, 2, ParamSpaceTag>;  // 1D -> 2D param
```

Progress:
✓ Updated Path3D to properly inherit from Manifold
✓ Converted to stack-based storage for derivatives
✓ Fixed inheritance hierarchy and virtual methods
✓ Basic tests passing

Remaining steps:
- Convert remaining path implementations (GeodesicCurve, PathSegment, TransitionPath)
- Ensure path length calculations work with new framework
- Add path-specific manifold tests

### 2. Collections and Connections
```cpp
template<int ParamDim, int TargetDim, typename SpaceTag>
class ManifoldCollection {
    // Will unify SurfaceCollection and path transition logic
};
```

Key steps:
- Create ManifoldCollection template
- Migrate existing collection functionality
- Implement unified connection system
- Add connection validation tests

## Implementation Plan

1. Path Implementation
   - Convert Path3D base class
   - Update path implementations
   - Verify path length calculations
   - Add path manifold tests

2. Collections
   - Create ManifoldCollection template
   - Migrate surface collection logic
   - Add path transition support
   - Implement connection validation
   - Add collection tests

## Benefits
1. Unified infrastructure for paths and surfaces
2. Better type safety through templates
3. Simplified maintenance with shared code
4. Cleaner architecture for future extensions

## Testing Strategy
- Convert path tests to use manifold interface
- Add specific tests for path length invariants
- Verify connection behavior in collection tests
- Ensure backward compatibility maintained

## Files to Update
- include/shap/path.hpp
- include/shap/path3d.hpp
- src/path.cpp
- src/path3d.cpp
- tests/path_tests.cpp
- tests/path_length_tests.cpp
- tests/manifold_tests.cpp

## New Files Needed
- include/shap/manifold_collection.hpp
- src/manifold_collection.cpp
- tests/manifold_collection_tests.cpp
