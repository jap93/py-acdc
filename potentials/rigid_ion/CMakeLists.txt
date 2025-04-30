cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(rigidion LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -O2 -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-Werror -Og")


find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(pybind11 REQUIRED)

set(SOURCES ${SOURCES} 
    src/Field.cpp
    src/TwoBody.cpp  
    src/NbrListPBC.cpp  
    src/Species.cpp
    src/Ewald.cpp
    src/Utility.cpp
)

add_library(rigidion MODULE "${SOURCES}")
target_link_libraries(rigidion PRIVATE pybind11::module Eigen3::Eigen )
set_target_properties(rigidion PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                     SUFFIX "${PYTHON_MODULE_EXTENSION}")
