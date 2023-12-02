cmake_minimum_required(VERSION 3.21)
project(SlauTest)

set(CMAKE_CXX_STANDARD 17)

set(SLAU_EXE_FILES
    GaussExe.cpp
    ZeidelExe.cpp
    LUdecompExe.cpp
    SimpleIterExe.cpp
    GradientsExe.cpp
)

foreach(EXE ${SLAU_EXE_FILES})
    string(REPLACE ".cpp" "" EXE_NAME ${EXE})
    add_executable(${EXE_NAME} ${EXE})
    target_link_libraries(${EXE_NAME} PUBLIC simfor)
endforeach()


#TODO
# find_package(OpenMP REQUIRED)
# if (OpenMP_CXX_FOUND)
# endif(OpenMP_CXX_FOUND)