#pragma once
#include "validation_config.hpp"
#include <array>
#include <cmath>
#include <stdexcept>

namespace shap {

// Tags for coordinate type
struct PointTag {};
struct VectorTag {};

// Tags for coordinate space
struct WorldSpaceTag {};
struct ParamSpaceTag {};

/**
 * Template class for N-dimensional coordinates.
 * 
 * @tparam N Dimensionality (2 or 3)
 * @tparam CoordTag PointTag or VectorTag
 * @tparam SpaceTag WorldSpaceTag or ParamSpaceTag
 */
template<int N, typename CoordTag, typename SpaceTag>
class Coord {
    static_assert(N >= 1 && N <= 3, "Only 1D, 2D, and 3D coordinates are supported");
    static_assert(std::is_same_v<CoordTag, PointTag> || std::is_same_v<CoordTag, VectorTag>,
                 "CoordTag must be either PointTag or VectorTag");
    static_assert(std::is_same_v<SpaceTag, WorldSpaceTag> || std::is_same_v<SpaceTag, ParamSpaceTag>,
                 "SpaceTag must be either WorldSpaceTag or ParamSpaceTag");

protected:
    std::array<double, N> coords_;

    // Copy constructor that allows conversion between Point and Vector types
    template<typename OtherTag>
    explicit Coord(const Coord<N, OtherTag, SpaceTag>& other) noexcept
        : coords_(other.coords_) {}

public:
    // Allow access to coords_ from other Coord instantiations
    template<int M, typename CT, typename ST>
    friend class Coord;

    using ThisType = Coord<N, CoordTag, SpaceTag>;
    using VectorType = Coord<N, VectorTag, SpaceTag>;  // Corresponding vector type

    /**
     * Default constructor - initializes all components to zero.
     */
    Coord() noexcept : coords_{} {}  // Zero-initialize array

    /**
     * Construct from individual components.
     */
    template<typename... Args>
    explicit Coord(Args... args) noexcept 
        : coords_{static_cast<double>(args)...} {
        static_assert(sizeof...(Args) == N, "Must provide exactly N components");
    }

    // Component access
    [[nodiscard]] double& operator[](int i) noexcept { return coords_[i]; }
    [[nodiscard]] double operator[](int i) const noexcept { return coords_[i]; }

    // World space accessors
    [[nodiscard]] double x() const noexcept requires std::is_same_v<SpaceTag, WorldSpaceTag> { return coords_[0]; }
    [[nodiscard]] double y() const noexcept requires std::is_same_v<SpaceTag, WorldSpaceTag> { return coords_[1]; }
    [[nodiscard]] double z() const noexcept requires (std::is_same_v<SpaceTag, WorldSpaceTag> && N == 3) { return coords_[2]; }

    // Parameter space accessors
    [[nodiscard]] double u() const noexcept requires std::is_same_v<SpaceTag, ParamSpaceTag> { return coords_[0]; }
    [[nodiscard]] double v() const noexcept requires std::is_same_v<SpaceTag, ParamSpaceTag> { return coords_[1]; }
    [[nodiscard]] double w() const noexcept requires (std::is_same_v<SpaceTag, ParamSpaceTag> && N == 3) { return coords_[2]; }

    // Conversion from 3D to 2D
    [[nodiscard]] Coord<2, CoordTag, SpaceTag> xy() const noexcept 
        requires (N == 3 && std::is_same_v<SpaceTag, WorldSpaceTag>) {
        return Coord<2, CoordTag, SpaceTag>(coords_[0], coords_[1]);
    }

    [[nodiscard]] Coord<2, CoordTag, SpaceTag> uv() const noexcept 
        requires (N == 3 && std::is_same_v<SpaceTag, ParamSpaceTag>) {
        return Coord<2, CoordTag, SpaceTag>(coords_[0], coords_[1]);
    }

    // Basic arithmetic for all coordinates
    ThisType operator+(const ThisType& other) const noexcept {
        ThisType result(*this);
        for (int i = 0; i < N; ++i) {
            result.coords_[i] += other.coords_[i];
        }
        return result;
    }

    VectorType operator-(const ThisType& other) const noexcept {
        VectorType result(*this);
        for (int i = 0; i < N; ++i) {
            result.coords_[i] -= other.coords_[i];
        }
        return result;
    }

    ThisType operator*(double scale) const noexcept {
        ThisType result(*this);
        for (int i = 0; i < N; ++i) {
            result.coords_[i] *= scale;
        }
        return result;
    }

    // Vector operations (available for vectors only)
    template<typename T = CoordTag>
    [[nodiscard]] double length_squared() const noexcept requires std::is_same_v<T, VectorTag> {
        double sum = 0.0;
        for (int i = 0; i < N; ++i) {
            sum += coords_[i] * coords_[i];
        }
        return sum;
    }

    template<typename T = CoordTag>
    [[nodiscard]] double length() const noexcept requires std::is_same_v<T, VectorTag> {
        return std::sqrt(length_squared());
    }

    template<typename T = CoordTag>
    [[nodiscard]] ThisType normalized() const requires std::is_same_v<T, VectorTag> {
        const double len = length();
        if (len < ValidationConfig::instance().vector_length_epsilon()) {
            throw std::invalid_argument("Cannot normalize zero-length vector");
        }
        ThisType result(*this);
        for (int i = 0; i < N; ++i) {
            result.coords_[i] /= len;
        }
        return result;
    }

    // Dot product (available for vectors only)
    template<typename T = CoordTag>
    [[nodiscard]] double dot(const ThisType& other) const noexcept requires std::is_same_v<T, VectorTag> {
        double sum = 0.0;
        for (int i = 0; i < N; ++i) {
            sum += coords_[i] * other.coords_[i];
        }
        return sum;
    }

    // Cross product (available for 3D vectors only)
    template<typename T = CoordTag>
    [[nodiscard]] ThisType crossed(const ThisType& other) const noexcept 
        requires (std::is_same_v<T, VectorTag> && N == 3) {
        return ThisType(
            coords_[1] * other.coords_[2] - coords_[2] * other.coords_[1],
            coords_[2] * other.coords_[0] - coords_[0] * other.coords_[2],
            coords_[0] * other.coords_[1] - coords_[1] * other.coords_[0]
        );
    }
};

// Point-specific operations
template<int N, typename SpaceTag>
[[nodiscard]] Coord<N, VectorTag, SpaceTag> operator-(
    const Coord<N, PointTag, SpaceTag>& a,
    const Coord<N, PointTag, SpaceTag>& b
) noexcept {
    Coord<N, VectorTag, SpaceTag> result(a);
    for (int i = 0; i < N; ++i) {
        result[i] -= b[i];
    }
    return result;
}

template<int N, typename SpaceTag>
[[nodiscard]] Coord<N, PointTag, SpaceTag> operator+(
    const Coord<N, PointTag, SpaceTag>& p,
    const Coord<N, VectorTag, SpaceTag>& v
) noexcept {
    Coord<N, PointTag, SpaceTag> result(p);
    for (int i = 0; i < N; ++i) {
        result[i] += v[i];
    }
    return result;
}

template<int N, typename SpaceTag>
[[nodiscard]] Coord<N, PointTag, SpaceTag> operator-(
    const Coord<N, PointTag, SpaceTag>& p,
    const Coord<N, VectorTag, SpaceTag>& v
) noexcept {
    Coord<N, PointTag, SpaceTag> result(p);
    for (int i = 0; i < N; ++i) {
        result[i] -= v[i];
    }
    return result;
}

// Vector-specific operations
template<int N, typename SpaceTag>
[[nodiscard]] Coord<N, VectorTag, SpaceTag> operator*(
    double scale,
    const Coord<N, VectorTag, SpaceTag>& v
) noexcept {
    return v * scale;
}

// Type aliases for common coordinate types
// 1D coordinates (for paths)
using ParamPoint1 = Coord<1, PointTag, ParamSpaceTag>;
using ParamVector1 = Coord<1, VectorTag, ParamSpaceTag>;

// 2D coordinates
using ParamPoint2 = Coord<2, PointTag, ParamSpaceTag>;
using ParamVector2 = Coord<2, VectorTag, ParamSpaceTag>;

// 3D coordinates
using WorldPoint3 = Coord<3, PointTag, WorldSpaceTag>;
using WorldVector3 = Coord<3, VectorTag, WorldSpaceTag>;
using ParamPoint3 = Coord<3, PointTag, ParamSpaceTag>;
using ParamVector3 = Coord<3, VectorTag, ParamSpaceTag>;

} // namespace shap
