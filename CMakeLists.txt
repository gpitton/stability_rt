project ("stability")
message(${PROJECT_SOURCE_DIR})
cmake_minimum_required (VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

#include_directories( INTERFACE /home/giuseppe/bin/Eigen/Dense)

add_executable(stability stability.cc
                         utils.cc
                         basis.cc
                         darcy.cc
                         )
