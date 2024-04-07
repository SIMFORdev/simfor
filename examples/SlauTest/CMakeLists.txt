cmake_minimum_required(VERSION 3.21)
project(SlauTest)

# find_package(MPI REQUIRED)
# find_package(Boost 1.84.0 REQUIRED)
add_compile_options(-O3)

set(CMAKE_CXX_STANDARD 17)

set(SLAU_EXE_FILES
    GaussExe.cpp
    GaussOmpExe.cpp
    GaussMpiExe.cpp
    ZeidelExe.cpp
    ZeidelOmpExe.cpp
    ZeidelMpiExe.cpp
    LUdecompExe.cpp
    LUdecompOmpExe.cpp
    SimpleIterExe.cpp
    SimpleIterOmpExe.cpp
    SimpleIterMpiExe.cpp
    GradientsExe.cpp
    GradientsOmpExe.cpp
    GradientsMpiExe.cpp
    TridiagonalExe.cpp
    TridiagonalOmpExe.cpp
    TridiagonalMpiExe.cpp
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