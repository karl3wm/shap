#include "coord.hpp"
#pragma once
#include <cstdint>

namespace shap {

/**
 * Represents indices for parameters in surface parameter space.
 * Used to identify which parameter (u or v) is being referenced.
 */
enum class ParamIndex : uint8_t {
    U = 0,  // First parameter coordinate
    V = 1   // Second parameter coordinate
};

// Arithmetic operators for parameter indices
inline int operator-(ParamIndex a, ParamIndex b) {
    return static_cast<int>(a) - static_cast<int>(b);
}

inline ParamIndex operator+(ParamIndex a, int b) {
    return static_cast<ParamIndex>(static_cast<int>(a) + b);
}

inline ParamIndex operator-(ParamIndex a, int b) {
    return static_cast<ParamIndex>(static_cast<int>(a) - b);
}

} // namespace shap
