cmake_minimum_required(VERSION 3.9)
project(CNSRV)

set(CMAKE_CXX_STANDARD 11)

add_executable(CNSRV
        Clustering.cpp
        Clustering.h
        main.cpp
        Orthology.cpp
        Orthology.h
        SpeciesNetwork.cpp
        SpeciesNetwork.h)
