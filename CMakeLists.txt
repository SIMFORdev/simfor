cmake_minimum_required(VERSION 3.21)
project(simfor
        VERSION 0.0.1
        DESCRIPTION "SIMFOR framework library"
        LANGUAGES ASM C CXX)

option(SIMFOR_BUILD_TESTS "Build SIMFOR tests" OFF)
option(SIMFOR_BUILD_EXAMPLES "Build SIMFOR examples" ON)
option(SIMFOR_BUILD_DOCS "Build SIMFOR documentation" OFF)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost 1.84.0 REQUIRED)

include_directories(${Boost_INCLUDE_DIR})

add_library(${PROJECT_NAME})

# target_link_libraries(${PROJECT_NAME} PUBLIC boost_mpi OpenMP::OpenMP_CXX MPI::MPI_CXX)

add_compile_options(-O3)
# target_compile_options(${PROJECT_NAME} PUBLIC ${OpenMP_CXX_FLAGS})
set(SIMFOR_INCLUDES
        # Public API includes
        include/simfor/SomeClass.hpp
        include/simfor/Gauss.hpp
        # include/simfor/GaussOmp.hpp
        # include/simfor/GaussMpi.hpp
        include/simfor/Zeidel.hpp
        # include/simfor/ZeidelOmp.hpp
        # include/simfor/ZeidelMpi.hpp
        include/simfor/LUdecomp.hpp        
        # include/simfor/LUdecompOmp.hpp        
        include/simfor/SimpleIter.hpp
        # include/simfor/SimpleIterOmp.hpp
        # include/simfor/SimpleIterMpi.hpp
        include/simfor/Gradients.hpp
        # include/simfor/GradientsOmp.hpp
        # include/simfor/GradientsMpi.hpp
        include/simfor/Tridiagonal.hpp
        # include/simfor/TridiagonalOmp.hpp
        # include/simfor/TridiagonalMpi.hpp
        # Private API includes
        include/simfor/internal/types.hpp
)

set(SIMFOR_SOURCES
        ${SIMFOR_INCLUDES}
        src/SomeClass.cpp
        src/Gauss.cpp
        # src/GaussOmp.cpp
        # src/GaussMpi.cpp
        src/Zeidel.cpp
        # src/ZeidelOmp.cpp
        # src/ZeidelMpi.cpp
        src/LUdecomp.cpp
        # src/LUdecompOmp.cpp
        src/SimpleIter.cpp
        # src/SimpleIterOmp.cpp
        # src/SimpleIterMpi.cpp
        src/Gradients.cpp
        # src/GradientsOmp.cpp
        # src/GradientsMpi.cpp
        src/Tridiagonal.cpp
        # src/TridiagonalOmp.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${SIMFOR_SOURCES})

target_sources(${PROJECT_NAME} PRIVATE ${SIMFOR_SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/${PROJECT_NAME}>
        $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
)

target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

if (SIMFOR_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif ()

if (SIMFOR_BUILD_EXAMPLES)
    include(cmake/utils.cmake)
    set(EXAMPLES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/examples")
    subdirlist(SUBDIRS ${EXAMPLES_DIR})
    foreach (subdir ${SUBDIRS})
        add_subdirectory("${EXAMPLES_DIR}/${subdir}")
    endforeach ()
endif ()

if (SIMFOR_BUILD_DOCS)
    include(cmake/doxygen.cmake)
endif ()