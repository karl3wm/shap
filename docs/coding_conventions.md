# Coding Conventions

## Core Principles

1. **Generalization Over Specialization**
   - Identify and lift specific solutions to their most general form
   - Create abstractions that capture underlying patterns
   - Prefer widely applicable solutions over context-specific ones

2. **Composition Over Complexity**
   - Build complex behavior from simple, well-defined components
   - Design interfaces that can be composed in predictable ways
   - Make components reusable across different contexts

3. **Strong Types Over Raw Values**
   - Use types to encode semantic meaning
   - Let the type system enforce invariants
   - Make invalid states unrepresentable

4. **Explicit Over Implicit**
   - Make relationships and transformations visible in the code
   - Document assumptions and invariants
   - Use clear naming to indicate purpose and context
   - Example: Prefix methods and variables with their space context (world_/parameter_)
   - Use shortened forms for common geometric concepts (e.g., 'pos' for 'position', 'norm' for 'normal')

5. **Brevity Scales with Usage**
   - Names should be shorter for more fundamental/frequently used constructs
   - Longer names are justified only for specialized/rarely used components
   - Reduce cognitive load for common operations

## Method Naming

1. **Method Names and Mutation**
   - When a method name could be interpreted as either a command or a getter:
     - Use imperative verb forms (e.g., normalize()) for methods that mutate the object
     - Use adjective/noun forms (e.g., normalized(), length()) for const methods that return values
   - This makes the mutation behavior immediately clear from the method name:
     - Commands to act -> mutation
     - Properties/states -> const getters

## Epsilon Values

1. **Validation Epsilons**: Used for validating geometric properties. Access via `ValidationConfig::instance()`.

2. **Algorithm Epsilons**: Used in specific algorithms as explicit parameters to maintain algorithm visibility and allow user control.

## Manifolds and Coordinate Spaces

The codebase uses a unified manifold-based approach for geometric objects:

1. **Manifold Template Parameters**
   - `ParamDim`: Dimension of parameter space (1 for paths, 2/3 for surfaces)
   - `TargetDim`: Dimension of target space (2 or 3 for world space)
   - `SpaceTag`: Indicates target space type (WorldSpaceTag or ParamSpaceTag)

2. **Parameter Space**
   - Domain: [0,1]^n where n is ParamDim
   - Type: `Coord<ParamDim, PointTag, ParamSpaceTag>`
   - Used for local coordinates on manifolds

3. **Target Space**
   - Domain: ℝ^n where n is TargetDim
   - Type: `Coord<TargetDim, PointTag, SpaceTag>`
   - Can be world space or parameter space of another manifold

4. **Common Manifold Types**
   - `WorldPath3D`: 1D → 3D world space paths
   - `Surface3D`: 2D → 3D world space surfaces
   - `ParamPath2D`: 1D → 2D parameter space paths
   - See manifold.hpp for complete list

## Error Handling

- Define errors in terms of violated invariants
- Provide clear context when invariants fail
- Use type system to prevent errors where possible

## Testing

- Test general properties that should hold universally
- Verify composition of operations maintains invariants
- Use shared utilities to express common test patterns
