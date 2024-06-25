cmake_minimum_required(VERSION 3.21)
project(MatrixTest)

set(CMAKE_CXX_STANDARD 17)

set(MATRIX_EXE_FILES
	classicEXE.cpp
	classic_ompEXE.cpp
	strassenEXE.cpp
	strassen_ompEXE.cpp
	strassen_recursiveEXE.cpp
	cblockEXE.cpp
	strassen_recursive_ompEXE.cpp
	multiplyVectorMatrixEXE.cpp
)

foreach(EXE ${MATRIX_EXE_FILES})
    string(REPLACE ".cpp" "" EXE_NAME ${EXE})
    add_executable(${EXE_NAME} ${EXE})
    target_link_libraries(${EXE_NAME} PUBLIC simfor)
endforeach()


#TODO
# find_package(OpenMP REQUIRED)
# if (OpenMP_CXX_FOUND)
# endif(OpenMP_CXX_FOUND)
