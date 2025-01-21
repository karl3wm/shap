cmake_minimum_required(VERSION 3.12)
project(shap VERSION 0.1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add include directory
include_directories(${PROJECT_SOURCE_DIR}/include)

# Create library target with implementation files
add_library(shap 
    src/surface.cpp
    src/path.cpp
    src/metric.cpp
)
target_include_directories(shap PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Build examples
add_executable(basic_demo examples/basic_demo.cpp)
target_link_libraries(basic_demo PRIVATE shap)

# Install rules
install(TARGETS shap
    EXPORT shapTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(DIRECTORY include/shap
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)

install(EXPORT shapTargets
    FILE shapTargets.cmake
    NAMESPACE shap::
    DESTINATION lib/cmake/shap
)
