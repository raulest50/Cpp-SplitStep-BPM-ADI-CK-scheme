cmake_minimum_required(VERSION 3.31)
project(CppSplitStepBPMADICKScheme)

set(CMAKE_CXX_STANDARD 20)

add_executable(CppSplitStepBPMADICKScheme main.cpp
        LaserBeam.cpp
        LaserBeam.h
        PrintHelpers.cpp
        PrintHelpers.h
        custom_solvers.cpp
        custom_solvers.h)

target_include_directories(CppSplitStepBPMADICKScheme
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)