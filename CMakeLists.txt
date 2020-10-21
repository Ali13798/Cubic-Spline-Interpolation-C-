cmake_minimum_required(VERSION 3.17)
project(Cubic_Spline_Interpolation)

set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

add_executable(Cubic_Spline_Interpolation main.cpp Matrix.cpp Matrix.h General_Include.h Invertible_Matrices.cpp Invertible_Matrices.h)