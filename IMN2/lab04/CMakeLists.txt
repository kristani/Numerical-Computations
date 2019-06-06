cmake_minimum_required(VERSION 3.12)
project(lab04)

set(CMAKE_CXX_STANDARD 14)
set ( CMAKE_CXX_FLAGS_RELEASE)
add_executable(lab04 main.cpp Node.h Element.cpp Element.h Equation.cpp Equation.h)