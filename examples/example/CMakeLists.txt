cmake_minimum_required(VERSION 3.23)
project(example)

set(CMAKE_CXX_STANDARD 17)

add_executable(${PROJECT_NAME} example.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC simfor)