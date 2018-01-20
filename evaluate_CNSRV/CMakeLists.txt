cmake_minimum_required(VERSION 3.5)
project(evaluate_CNSRV)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -O3")

set(SOURCE_FILES main.cpp)
add_executable(evaluate_CNSRV ${SOURCE_FILES})

