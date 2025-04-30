cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(metal LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -O2 -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-Werror -Og")


find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(pybind11 REQUIRED)

set(SOURCES ${SOURCES} 
    src/Field.cpp
    src/ManyBody.cpp  
    src/NbrListPBC.cpp  
    src/Species.cpp
    src/MetalPotential.cpp
)

add_library(metal MODULE "${SOURCES}")
target_link_libraries(metal PRIVATE pybind11::module Eigen3::Eigen )
set_target_properties(metal PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                     SUFFIX "${PYTHON_MODULE_EXTENSION}")
