#include "coord.hpp"
#ifndef SHAP_VALIDATION_CONFIG_HPP
#define SHAP_VALIDATION_CONFIG_HPP

namespace shap {

/**
 * Global configuration singleton for validation-related epsilon values.
 * These epsilons are used for validating geometric properties like
 * vector lengths and parallel vector checks.
 */
class ValidationConfig {
public:
    static ValidationConfig& instance() {
        static ValidationConfig instance;
        return instance;
    }

    // Epsilon for checking if a vector has zero length
    double vector_length_epsilon() const { return vector_length_epsilon_; }
    void set_vector_length_epsilon(double epsilon) { vector_length_epsilon_ = epsilon; }

    // Epsilon for checking if vectors are parallel
    double vector_parallel_epsilon() const { return vector_parallel_epsilon_; }
    void set_vector_parallel_epsilon(double epsilon) { vector_parallel_epsilon_ = epsilon; }

    // Epsilon for parameter bound checks
    double parameter_bound_epsilon() const { return parameter_bound_epsilon_; }
    void set_parameter_bound_epsilon(double epsilon) { parameter_bound_epsilon_ = epsilon; }

private:
    ValidationConfig() 
        : vector_length_epsilon_(1e-10)
        , vector_parallel_epsilon_(1e-10)
    {}

    // Delete copy/move operations to ensure singleton
    ValidationConfig(const ValidationConfig&) = delete;
    ValidationConfig& operator=(const ValidationConfig&) = delete;
    ValidationConfig(ValidationConfig&&) = delete;
    ValidationConfig& operator=(ValidationConfig&&) = delete;

    double vector_length_epsilon_;
    double vector_parallel_epsilon_;
    double parameter_bound_epsilon_ = 1e-10;
};

} // namespace shap

#endif // SHAP_VALIDATION_CONFIG_HPP
