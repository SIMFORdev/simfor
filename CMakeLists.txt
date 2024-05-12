cmake_minimum_required(VERSION 3.10)
project(Plotter)

# Source files
set(SOURCES
    src/plotter.cpp
    examples/test_2d.cpp
    examples/test_3d.cpp
    examples/test_points.cpp
)

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# Include and link directories
include_directories(${CMAKE_SOURCE_DIR}/src)
link_directories(${CMAKE_BINARY_DIR}/lib)

# Library
add_library(plotter STATIC src/plotter.cpp)

# Executables
add_executable(test_2d examples/test_2d.cpp)
add_executable(test_3d examples/test_3d.cpp)
add_executable(test_points examples/test_points.cpp)

# Link libraries to executables
target_link_libraries(test_2d plotter)
target_link_libraries(test_2d GL glut)

target_link_libraries(test_3d plotter)
target_link_libraries(test_3d GL glut)

target_link_libraries(test_points plotter)
target_link_libraries(test_points GL glut)
