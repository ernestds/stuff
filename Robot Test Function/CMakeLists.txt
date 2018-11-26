cmake_minimum_required(VERSION 3.0)

project(Test)

SET(GCC_COVERAGE_COMPILE_FLAGS "-O2 -pthread -fopenmp -fopenmp-simd") 
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB program_sources ./sourcefiles/*.cpp ./sourcefiles/*/*.cpp ./sourcefiles/*/*/*.cpp ./sourcefiles/*/*/*/*.cpp ./sourcefiles/*.c ./sourcefiles/*/*.c ./sourcefiles/*/*/*.c ./sourcefiles/*/*/*/*.c)
SET( program ${program_sources})

find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories( ${EIGEN3_INCLUDE_DIRS} )

ADD_EXECUTABLE(program.o ${program})