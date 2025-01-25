#include "shap/path3d.hpp"
#include <stdexcept>

namespace shap {

WorldVector3
Path3D::normal(const ParamPoint1& param) const {
    return normal_(param);
}

WorldVector3
Path3D::binormal(const ParamPoint1& param) const {
    return tangent_(param).crossed(normal_(param));
}

} // namespace shap
