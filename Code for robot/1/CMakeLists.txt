cmake_minimum_required(VERSION 3.0)

project(Test)

list(APPEND CMAKE_MODULE_PATH /home/franka/Projects/libfranka/cmake)

SET(GCC_COVERAGE_COMPILE_FLAGS "-O2 -pthread") 

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Franka REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Poco REQUIRED COMPONENTS Foundation)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

file(GLOB program_sources ./sourcefiles/*.cpp ./sourcefiles/*/*.cpp ./sourcefiles/*/*/*.cpp ./sourcefiles/*/*/*/*.cpp ./sourcefiles/*.c ./sourcefiles/*/*.c ./sourcefiles/*/*/*.c ./sourcefiles/*/*/*/*.c)
SET( program ${program_sources})

add_library(examples_common STATIC
	examples_common.cpp
)

target_link_libraries(examples_common PUBLIC Franka::Franka Eigen3::Eigen3)

find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories( ${EIGEN3_INCLUDE_DIRS} )


add_executable(program.o ${program})
target_link_libraries(program.o examples_common Franka::Franka Eigen3::Eigen3 Threads::Threads)
