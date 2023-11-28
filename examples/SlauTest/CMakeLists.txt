cmake_minimum_required(VERSION 3.21)
project(SlauTest)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

if (OpenMP_CXX_FOUND)
    add_executable(GaussOmp GaussOmpExe.cpp)
    target_link_libraries(GaussOmp PUBLIC simfor OpenMP::OpenMP_CXX)
endif(OpenMP_CXX_FOUND)
