cmake_minimum_required(VERSION 3.23)
project(simfor
        VERSION 0.0.1
        DESCRIPTION "SIMFOR framework library"
        LANGUAGES ASM C CXX)

option(SIMFOR_BUILD_TESTS "Build SIMFOR tests" OFF)
option(SIMFOR_BUILD_EXAMPLES "Build SIMFOR examples" ON)
option(SIMFOR_BUILD_DOCS "Build SIMFOR documentation" OFF)

add_library(${PROJECT_NAME})

add_compile_options(-O3)

set(SIMFOR_INCLUDES
        # Public API includes
        include/simfor/SomeClass.hpp
        # Private API includes
        include/simfor/internal/types.hpp
)

set(SIMFOR_SOURCES
        ${SIMFOR_INCLUDES}
        src/SomeClass.cpp
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