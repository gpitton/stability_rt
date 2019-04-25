project ("stability")
cmake_minimum_required (VERSION 2.6)

#set (CMAKE_CXX_FLAGS "-DEIGEN_STACK_ALLOCATION_LIMIT=0 -Wall -pedantic")
#set (CMAKE_CXX_FLAGS "-DEIGEN_STACK_ALLOCATION_LIMIT=0 -O3")
#set(CMAKE_CXX_FLAGS "-DEIGEN_STACK_ALLOCATION_LIMIT=0 -fopenmp -Wall -pedantic")
set(CMAKE_CXX_FLAGS "-DEIGEN_STACK_ALLOCATION_LIMIT=0 -fopenmp -O3 -std=gnu++17")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#set(Spectra_INCLUDE_DIR /home/giuseppe/bin/Spectra-0.8.0)
#set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
#set(Boost_USE_STATIC_RUNTIME OFF) 

find_package(Boost 1.45.0 REQUIRED)
find_package(Eigen3 REQUIRED)
#find_package(Spectra REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Spectra_INCLUDE_DIR})

add_executable(stability stability.cc)
add_executable(two_parameters two_parameters.cc)
