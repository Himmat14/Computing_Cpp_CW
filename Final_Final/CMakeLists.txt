cmake_minimum_required(VERSION 3.28)
project(Final_Final)

set(CMAKE_CXX_STANDARD 17)

add_executable(Final_Final main.cpp
        cmake-build-debug/Simulation.cpp
        cmake-build-debug/Simulation.h
        cmake-build-debug/EulerSimulation.cpp
        cmake-build-debug/EulerSimulation.h
        cmake-build-debug/RK4Simulation.cpp
        cmake-build-debug/RK4Simulation.h
        cmake-build-debug/Vessel.cpp
        cmake-build-debug/Vessel.h)
