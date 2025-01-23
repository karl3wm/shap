# Coding Conventions

## Epsilon Values

There are two categories of epsilon values used in the codebase:

1. **Validation Epsilons**: Used for validating geometric properties like vector lengths and parallel vector checks. These are managed through the `ValidationConfig` singleton to ensure consistent validation across the codebase. Access these via `ValidationConfig::instance()`. Current validation epsilons include:
   - vector_length_epsilon: For checking if a vector has zero length
   - vector_parallel_epsilon: For checking if vectors are parallel

2. **Algorithm Epsilons**: Used in specific algorithms (like projection, intersection calculations, or parameter bound checks). These should remain as explicit parameters to maintain algorithm visibility and allow user control. Examples include:
   - Degenerate case handling in coordinate transformations (e.g., world_to_parameter_space_with_epsilon)
   - Parameter bound checks in path solving
   - Intersection tolerance checks

This separation ensures consistent validation behavior while preserving explicit control over algorithm-specific numerical tolerances. Validation epsilons are hidden implementation details, while algorithm epsilons are first-class parameters that users can tune based on their needs.


## Core Principles

1. **Generalization Over Specialization**
   - Identify and lift specific solutions to their most general form
   - Create abstractions that capture underlying patterns
   - Prefer widely applicable solutions over context-specific ones
   - Example: Instead of specialized test utilities per file, create general-purpose testing tools

2. **Composition Over Complexity**
   - Build complex behavior from simple, well-defined components
   - Design interfaces that can be composed in predictable ways
   - Make components reusable across different contexts
   - Example: Surface transformations compose parameter space and world space operations

3. **Strong Types Over Raw Values**
   - Use types to encode semantic meaning
   - Let the type system enforce invariants
   - Make invalid states unrepresentable
   - Example: `ParameterCoordinate` enforces [0,1] bounds instead of raw doubles

4. **Explicit Over Implicit**
   - Make relationships and transformations visible in the code
   - Document assumptions and invariants
   - Use clear naming to indicate purpose and context
   - Example: Prefix methods and variables with their space context (world_/parameter_)
   - Example: `world_u` for a vector in world space corresponding to parameter u direction
   - Remove redundant suffixes when the space context makes the type clear
   - Use shortened forms for common geometric concepts (e.g., 'pos' for 'position', 'norm' for 'normal')

5. **Brevity Scales with Usage**
   - Names should be shorter for more fundamental/frequently used constructs
   - Longer names are justified only for specialized/rarely used components
   - Reduce cognitive load for common operations
   - Example: `PlanarPatch` vs `TemporaryConfigurationStorage`

## Coordinate Spaces

The codebase operates on two fundamental spaces that exemplify our core principles:

1. **Parameter Space**
   - Domain: [0,1] × [0,1]
   - Strong types: `ParameterCoordinate`, `ParameterVelocity`
   - Enforced invariants: All values validated to be in [0,1]

2. **World Space**
   - Domain: ℝ³
   - Strong types: `Point`, `Vector`
   - Natural operations: Physical positions and directions

### Space Transformations

Space transformations demonstrate composition of our core principles:

1. **Type Safety**
   ```cpp
   // Types encode meaning and prevent mixing spaces
   ParameterCoordinate param(0.5, 0.5);
   Point world_pos = surface->evaluate(param).world_position();
   ```

2. **Scale Factors**
   ```cpp
   // Compose simple operations for complex transformations
   auto [du_scale, dv_scale] = surface->get_scale_factors(param);
   auto world_length = parameter_length * du_scale;
   ```

3. **Metric Tensor**
   ```cpp
   // General solution for space relationships
   auto param_vel = world_to_parameter_velocity(
       world_direction, world_du, world_dv);
   ```

## Implementation Guidelines

1. **Code Transformations**
   - Prefer using `sed` for systematic code changes
   - Create reusable transformation scripts for common patterns
   - Document transformations for future reference
   - Example: Use `sed` to update type names across files
   ```bash
   # Update Point to WorldPoint3 across files
   sed -i 's/\bPoint\b/WorldPoint3/g' include/shap/*.hpp src/*.cpp
   
   # Remove circular includes (e.g., a file including itself)
   sed -i '1{/^#include.*self\.hpp/d}' path/to/self.hpp
   
   # Update member variable names (e.g., removing underscores)
   sed -i 's/coords_/coords/g' path/to/file.hpp
   ```

   Always prefer using sed for systematic changes across files. This makes changes:
   - More consistent
   - Easier to review
   - Quicker to execute
   - Simpler to revert if needed

2. **Code Organization**
   - Group related functionality into composable units
   - Share common utilities across the codebase
   - Place general-purpose tools at appropriate scope levels
   - Example: Test utilities in shared namespace vs. file-specific helpers

2. **Documentation**
   - Document the general pattern, not just the specific use
   - Explain relationships between components
   - Show how specific cases follow from general principles
   - Example: Metric tensor documentation explains general space relationships

3. **Error Handling**
   - Define errors in terms of violated invariants
   - Provide clear context when invariants fail
   - Use type system to prevent errors where possible
   - Example: Parameter space bounds checking through types

4. **Testing**
   - Test general properties that should hold universally
   - Verify composition of operations maintains invariants
   - Use shared utilities to express common test patterns
   - Example: Path length preservation tests verify general metric properties

## Common Patterns

These patterns demonstrate how specific implementations follow from our core principles:

1. **Surface Point Creation**
```cpp
// Compose parameter and world space properties
SurfacePoint point(
    surface,
    ParameterCoordinate(0.5, 0.5),  // Strong type for parameters
    Point(1, 1, 0),                 // Strong type for position
    Vector(0, 0, 1),                // Strong type for direction
    Vector(2, 0, 0),                // Explicit scale factors
    Vector(0, 2, 0)
);
```

2. **Path Creation**
```cpp
// General pattern for creating paths on any surface
auto path = surface->create_path(
    start_point,
    world_direction,    // Explicit space context
    world_length       // Clear units
);
```

3. **Test Organization**
```cpp
// General utilities in shared namespace
namespace shap::test {
    // Common patterns lifted to general tools
    [[nodiscard]] constexpr bool approx_equal(double a, double b,
        double epsilon = EPSILON) noexcept;
}
```

## Best Practices

1. **Follow General Patterns**
   - Look for opportunities to generalize specific solutions
   - Reuse existing abstractions before creating new ones
   - Compose simple tools rather than building complex ones

2. **Maintain Invariants**
   - Use types to encode invariants where possible
   - Document invariants that cross component boundaries
   - Test that compositions preserve invariants

3. **Clear Communication**
   - Name things according to their general purpose
   - Document the general pattern being implemented
   - Show how specific uses follow from general principles
