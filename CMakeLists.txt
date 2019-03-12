# 
# Specify minimum CMake version and project name
#
cmake_minimum_required(VERSION 2.8)
project(project2)

#
# CMake setup
#
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set (CMAKE_VERBOSE_MAKEFILE 0) # 1 should be used for debugging
set (CMAKE_SUPPRESS_REGENERATION TRUE) # Suppress ZERO_CHECK
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O2 -lpthread -Wall -Wextra -Wpedantic -std=c++11")

find_package(MPI)

#
# Set include paths
#
include_directories (src/ ${MPI_INCLUDE_PATH})

#
# Add files
#
file (GLOB_RECURSE PROJECT_HEADERS src/*.hpp
								   src/*.h)
file (GLOB_RECURSE PROJECT_SOURCES src/*.cpp)
file (GLOB 		   PROJECT_CONFIGS CMakeLists.txt
								   README.rst
								   .gitignore
								   .gitmodules)

#
# Organizing files
#
source_group ("headers" FILES ${PROJECT_HEADERS})
source_group ("sources" FILES ${PROJECT_SOURCES})

#
# Set executable and target link libraries
#
add_definitions (-DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")
add_executable(${PROJECT_NAME} ${PROJECT_SOURCES} ${PROJECT_HEADERS} ${PROJECT_CONFIGS})
target_link_libraries (${PROJECT_NAME} ${MPI_C_LIBRARIES})
set_target_properties (${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${PROJECT_NAME})
