# Test Failure Investigation Documentation

This directory contains documentation for ongoing test failure investigations. Each document should follow these guidelines to ensure issues can be tracked and resolved effectively.

## Document Structure

Each test failure investigation document should include:

1. **Header Metadata**
   ```markdown
   ---
   status: [OPEN|RESOLVED|BLOCKED]
   test_file: path/to/test_file.cpp
   test_name: "Name of the failing test"
   first_seen: YYYY-MM-DD
   last_updated: YYYY-MM-DD
   resolution: [If resolved, brief description of fix]
   related_issues: [Links to related issues/PRs]
   ---
   ```

2. **Failure Description**
   - Exact test output showing expected vs actual values
   - Environment/configuration where failure occurs
   - Steps to reproduce
   - Any patterns in when the failure occurs

3. **Investigation Progress**
   - Current understanding of the issue
   - Hypotheses tested and results
   - Mathematical relationships being violated
   - Code areas identified as relevant
   - Debugging insights gained

4. **Resolution Status**
   - For OPEN: Current blockers and next steps
   - For RESOLVED: 
     - Root cause identified
     - Fix implemented
     - Verification steps performed
     - New tests added to prevent regression
   - For BLOCKED:
     - Dependencies or prerequisites needed
     - Impact assessment of the blockage

## Diagnostic Methodology

1. **Add Diagnostic Output**
   - Before attempting fixes, add detailed diagnostic output at key points
   - Log intermediate calculation steps and values
   - For mathematical operations (like metric tensors):
     * Log input values
     * Show step-by-step calculations
     * Compare with expected theoretical values
   - Use consistent formatting for easy comparison

2. **Analyze Output**
   - Compare diagnostic values against mathematical theory
   - Look for patterns in where values diverge
   - Verify each transformation step independently
   - Document findings in the investigation document

3. **Validate Understanding**
   - Create small, focused test cases that isolate behavior
   - Add assertions to verify intermediate calculations
   - Compare results with manual calculations
   - Document the validation process

4. **Only Then Fix**
   - Once the exact nature of the error is understood
   - When you can predict the impact of the fix
   - With clear evidence supporting the change
   - While preserving the diagnostic output

## Usage Guidelines

1. Create a new document when starting investigation of a test failure
2. Update the status and last_updated fields as investigation progresses
3. Document all significant findings and attempted solutions
4. Reference related code, PRs, and other documentation
5. When resolved, ensure the resolution is clearly documented for future reference
6. Follow the diagnostic methodology before attempting fixes

## Current Investigations

### Path Length Investigation

A series of test failures related to path length preservation and parameter space transformations:

1. **Path Length Validation** (path_length_validation.md)
   - Status: OPEN
   - Comprehensive analysis of path length calculations
   - Documents current test status:
     * Parameter/World mappings PASS
     * Path length preservation FAILS
   - Metric tensor computation verified correct
   - Issue isolated to path evaluation

2. **Length Scaling** (length_scaling.md)
   - Test: path_tests.cpp:test_cube_face_paths()
   - Status: OPEN
   - Issue: Path endpoints not matching expected positions
   - Key Finding: End points are off by 1.5x expected distance
   - Related: path_length_validation.md

3. **Space Transformations** (space_transformations.md)
   - Test: space_transformation_tests.cpp:test_space_transformations()
   - Status: OPEN
   - Results:
     * Basic coordinate mappings work correctly
     * Path distances show 1.5x error
   - Confirms issue is in path evaluation

4. **Parameter Space Analysis** (parameter_space.md, parameter_space_validation.md)
   - Status: OPEN
   - Focus: Parameter space relationships and validation
   - Dependencies: Requires resolution of path length issues

### Test Analysis Structure

The investigation follows this progression:

1. Coordinate Transformations (Verified Working)
   - Parameter to world space mapping
   - World to parameter space mapping
   - Basic geometric properties

2. Path Creation (Under Investigation)
   - Direction projection
   - Length scaling
   - Parameter velocity computation
   - Metric tensor application

3. Path Evaluation (Current Focus)
   - Distance preservation
   - Velocity consistency
   - Scale factor handling

Key Findings:
1. Basic geometry and transformations work correctly
2. Metric tensor computation is verified accurate
3. Issue appears in path evaluation where:
   - Parameter space calculations look correct
   - But world space distances are 1.5x too large
   - Suggesting a problem in distance accumulation
