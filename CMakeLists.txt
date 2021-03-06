project(poisson)

cmake_minimum_required(VERSION 2.6)

enable_language(CXX)
enable_language(Fortran)

# CMake Options to enable/disable MPI and OpenMP
option(ENABLE_OPENMP "Enable OpenMP support?" ON)
option(ENABLE_MPI    "Enable MPI support?" ON)

if(ENABLE_MPI)
  find_package(MPI)
endif()

if(MPI_FOUND)
  add_definitions(-DHAVE_MPI=1)
  set(INCLUDES ${INCLUDES} ${MPI_INCLUDE_PATH})
  set(DEPLIBS_CXX ${DEPLIBS_CXX} ${MPI_C_LIBRARIES})
  set(DEPLIBS_F ${DEPLIBS_F} ${MPI_Fortran_LIBRARIES})
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_definitions(-DHAVE_OPENMP=1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++11 -Wall -Wextra -Wpedantic")
  endif()
endif()

# Send the include directories to the compiler
include_directories(${INCLUDES})

# Add a library with common sources
add_library(common src/fst.f)

# The C++ based poisson application
add_executable(poisson src/poisson.cpp)
add_executable(parallelPoisson src/parallel-poisson.cpp)

target_link_libraries(poisson common ${DEPLIBS_CXX})
target_link_libraries(parallelPoisson common ${DEPLIBS_CXX})

set_target_properties(poisson PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(parallelPoisson PROPERTIES LINKER_LANGUAGE CXX)

## C, NOT USED ##




# The C based poisson application
# add_executable(poisson src/poisson.c)

# Link to libraries
#target_link_libraries(poisson common ${DEPLIBS_C})

# We want to link this as a C program
# set_target_properties(poisson PROPERTIES LINKER_LANGUAGE C)


## FORTRAN, WHICH WILL NOT BE USED ##

# The Fortran based poisson application
# add_executable(poisson-f src/poisson.f90)

# Link to libraries
# target_link_libraries(poisson-f common ${DEPLIBS_F})
