# Parameter Space Conventions

## Core Concepts

Every surface in the library uses a standardized parameter space [0,1]×[0,1]:
- Parameters are always normalized to the [0,1] range
- (u,v) coordinates uniquely identify points on the surface
- Surface edges correspond to parameter bounds (0 or 1)

## Surface Types

The library supports different types of surfaces with specific properties:

### Developable Surfaces (e.g., cube faces)
- Can be flattened without distortion
- Geodesics are straight lines when flattened
- Linear interpolation in parameter space gives exact paths
- Example: SquareFace maps (u,v) linearly to 3D space

### Smooth Surfaces (e.g., sphere)
- Have continuous derivatives
- Geodesics follow differential equations
- Paths require numerical integration
- Example: SphereSurface maps (u,v) to spherical coordinates

## Edge Connections

Surfaces can be connected along their edges:
- Each edge is identified by fixing one parameter (u or v) at 0 or 1
- The other parameter varies along the edge
- Connections specify how parameters align between surfaces
- Orientation (+1/-1) handles parameter reversal between surfaces

Example: Cube Face Connections
```
Front face (u=1) → Right face (u=0):
- Front face: y=1, x varies with u, z varies with v
- Right face: x=1, y varies with u, z varies with v
- Parameters increase in same direction (orientation = +1)
```

## Path Creation

Creating paths on surfaces follows these rules:
1. Start with a point (u,v) on the surface
2. Project desired direction onto surface
3. Convert world-space length to parameter space
4. Sample points along path maintaining surface constraints

## Path Length Handling

Path length must be carefully managed between world space and parameter space:

### World Space Length
- Specified in input as physical distance
- Measured along the actual 3D curve
- Independent of surface parameterization
- Example: Distance of 1.0 means move 1 unit in 3D space

### Parameter Space Length
- Must be scaled based on surface metrics
- For SquareFace:
  * du vector defines x-scale (e.g., 2 units wide)
  * dv vector defines z-scale (e.g., 2 units tall)
  * Parameter change of 1.0 maps to full surface width/height
  * Must scale parameter steps by surface dimensions

### Length Conversion
1. Project world direction onto surface
2. Compute parameter space direction
3. Scale parameter step by surface metric
4. Ensure total path length matches requested length

For developable surfaces like SquareFace:
- Linear interpolation in parameter space
- Gives exact straight lines in world space
- Must account for surface scale factors