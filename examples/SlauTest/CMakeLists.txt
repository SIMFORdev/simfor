cmake_minimum_required(VERSION 3.21)
project(SlauTest)

set(CMAKE_CXX_STANDARD 17)

set(SLAU_EXE_FILES
    GaussExe.cpp
    # GaussOmpExe.cpp
    # GaussMpiExe.cpp
    ZeidelExe.cpp
    # ZeidelOmpExe.cpp
    # ZeidelMpiExe.cpp
    LUdecompExe.cpp
    # LUdecompOmpExe.cpp
    SimpleIterExe.cpp
    # SimpleIterOmpExe.cpp
    # SimpleIterMpiExe.cpp
    GradientsExe.cpp
    # GradientsOmpExe.cpp
    # GradientsMpiExe.cpp
    TridiagonalExe.cpp
    # TridiagonalOmpExe.cpp
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