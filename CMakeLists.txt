cmake_minimum_required(VERSION 3.15)
project(shap VERSION 0.1.0 LANGUAGES CXX)

# Require C++20 for constexpr features and concepts
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-fconcepts)

# Add the include directory
include_directories(include)

# Add the example executable
add_executable(basic_demo examples/basic_demo.cpp)

# Add compile definitions for M_PI if needed
target_compile_definitions(basic_demo PRIVATE
    _USE_MATH_DEFINES
)
