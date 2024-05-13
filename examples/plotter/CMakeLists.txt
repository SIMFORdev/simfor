cmake_minimum_required(VERSION 3.21)
project(plotter)

#find_package(GLUT REQUIRED)
#find_package(OpenGL REQUIRED)

set(CMAKE_CXX_STANDARD 17)

set(PLOTTER_FILES
    test_2d.cpp
    test_3d.cpp
    test_points.cpp
)

foreach(EXE ${PLOTTER_FILES})
    string(REPLACE ".cpp" "" EXE_NAME ${EXE})
    add_executable(${EXE_NAME} ${EXE})
    target_link_libraries(${EXE_NAME} PUBLIC simfor)
endforeach()