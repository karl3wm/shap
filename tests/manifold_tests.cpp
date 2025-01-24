#include "shap/coord.hpp"
#include "shap/surface3d.hpp"
#include "shap/path3d.hpp"
#include "shap/surfaces/flat_patch.hpp"
#include "shap/geometric_point.hpp"
#include "shap/path.hpp"
#include "test_utils.hpp"
#include <cassert>
#include <iostream>

namespace shap::test {

void test_basic_manifold_operations() {
    // Create a flat patch as a simple 2D manifold
    auto face = surfaces::create_flat_patch(
        WorldPoint3(0, 0, 0),
        WorldVector3(1, 0, 0),
        WorldVector3(0, 1, 0)
    );

    // Test basic evaluation at center point
    const auto center_param = ParamPoint2(0.5, 0.5);
    const auto center_geom = face->evaluate(center_param);
    
    // Check position
    assert(approx_equal(center_geom.world_pos(), WorldPoint3(0.5, 0.5, 0)));
    
    // Check derivatives
    const auto derivs = face->derivatives(center_param);
    assert(derivs.size() == 2);  // Should have two derivatives for 2D manifold
    
    // First derivative should be (1,0,0) - along u direction
    assert(approx_equal(derivs[0], WorldVector3(1.0, 0.0, 0.0)));
    
    // Second derivative should be (0,1,0) - along v direction
    assert(approx_equal(derivs[1], WorldVector3(0.0, 1.0, 0.0)));
    
    std::cout << "Basic manifold operations test passed!\n";
}

} // namespace shap::test

int main() {
    try {
        std::cout << "Testing basic manifold operations...\n";
        shap::test::test_basic_manifold_operations();
        
        std::cout << "\nAll manifold interface tests passed!\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << "\n";
        return 1;
    }
}
