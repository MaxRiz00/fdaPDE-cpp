cmake_minimum_required(VERSION 3.14)

# Project name
project(SplineProject)

# Set C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable warnings and optimizations
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#set(CMAKE_WARNINGS_FLAGS "-Wall ")
set(CMAKE_CXX_FLAGS "-O2 -march=native -Wno-psabi")

# Find Eigen3
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

# Include directories
#include_directories(../fdaPDE) 
include_directories(../../fdaPDE/)
include_directories(../../fdaPDE/core/)

# Add the executable
add_executable(prova europe_solver.cpp)

# Link Eigen3 library
target_link_libraries(prova Eigen3::Eigen)
