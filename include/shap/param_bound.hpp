#include "coord.hpp"
#pragma once
#include <cstdint>

namespace shap {

/**
 * Represents bounds (0 or 1) in surface parameter space.
 * Used to identify which boundary of a parameter range is being referenced.
 */
enum class ParamBound : uint8_t {
    Lower = 0,  // Parameter value 0
    Upper = 1   // Parameter value 1
};

// Arithmetic operators for parameter bounds
inline int operator-(ParamBound a, ParamBound b) {
    return static_cast<int>(a) - static_cast<int>(b);
}

inline ParamBound operator+(ParamBound a, int b) {
    return static_cast<ParamBound>(static_cast<int>(a) + b);
}

inline ParamBound operator-(ParamBound a, int b) {
    return static_cast<ParamBound>(static_cast<int>(a) - b);
}

} // namespace shap
