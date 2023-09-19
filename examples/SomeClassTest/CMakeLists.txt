cmake_minimum_required(VERSION 3.21)
project(SomeClassTest)

set(CMAKE_CXX_STANDARD 17)

add_executable(${PROJECT_NAME} SomeClassTest.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC simfor)