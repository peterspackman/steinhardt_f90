cmake_minimum_required(VERSION 3.10)
project(steinhardt Fortran)

# Set the Fortran compiler and flags
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_Fortran_FLAGS "-Wall -Wextra -Wpedantic -std=f2008")

# Add the source files
set(SOURCES
    legendre.f90
    spherical_harmonics.f90
    steinhardt.f90
    wigner3j.f90
    main.f90
)

# Add the executable
add_executable(steinhardt ${SOURCES})
