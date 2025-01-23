#include "shap/coord.hpp"
#pragma once

#include <string_view>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace shap::test {

// Constants for numerical testing
constexpr double EPSILON = 1e-10;

// Helper to check if two doubles are approximately equal
[[nodiscard]] constexpr bool approx_equal(double a, double b, double epsilon = EPSILON) noexcept {
    return std::abs(a - b) <= epsilon;
}

// Helper to check if two points are approximately equal
[[nodiscard]] inline bool approx_equal(const WorldPoint3& a, const WorldPoint3& b, double epsilon = EPSILON) noexcept {
    const bool result = approx_equal(a.x(), b.x(), epsilon) &&
                       approx_equal(a.y(), b.y(), epsilon) &&
                       approx_equal(a.z(), b.z(), epsilon);
    
    if (!result) {
        std::cout << "Point comparison failed:\n"
                  << "  Expected: (" << a.x() << ", " << a.y() << ", " << a.z() << ")\n"
                  << "  Actual:   (" << b.x() << ", " << b.y() << ", " << b.z() << ")\n"
                  << "  Diff:     (" 
                  << std::abs(a.x() - b.x()) << ", "
                  << std::abs(a.y() - b.y()) << ", "
                  << std::abs(a.z() - b.z()) << ")\n"
                  << "  Epsilon:  " << epsilon << "\n";
    }
    return result;
}

// Print a point for debugging
inline void print_point(std::string_view label, const WorldPoint3& p) {
    std::cout << label << ": ("
              << std::fixed << std::setprecision(6)
              << p.x() << ", " << p.y() << ", " << p.z() << ")\n";
}

} // namespace shap::test
