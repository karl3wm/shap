#include "shap/path3d.hpp"
#include <stdexcept>

namespace shap {

namespace {
    void validate_parameters(const ParamPoint1& param) {
        if (param[0] < 0.0 || param[0] > 1.0) {
            throw std::invalid_argument("Path parameter must be in [0,1]");
        }
    }
}

WorldVector3
Path3D::normal(const ParameterPoint& param) const {
    const double t = get_param_value(param);
    return normal_(t);
}

WorldVector3
Path3D::binormal(const ParameterPoint& param) const {
    const double t = get_param_value(param);
    return tangent_(t).crossed(normal_(t));
}

double
Path3D::get_param_value(const ParameterPoint& param) {
    validate_parameters(param);
    return param[0];
}

} // namespace shap
