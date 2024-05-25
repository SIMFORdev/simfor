cmake_minimum_required(VERSION 3.8...3.21)
project(simfor
        VERSION 0.0.1
        DESCRIPTION "SIMFOR framework library"
        LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(utils)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

find_package(Boost 1.84.0 COMPONENTS mpi REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

option(SIMFOR_BUILD_TESTS "Build SIMFOR tests" OFF)
option(SIMFOR_BUILD_EXAMPLES "Build SIMFOR examples" ON)
option(SIMFOR_BUILD_DOCS "Build SIMFOR documentation" OFF)
option(SIMFOR_INSTALL "Generate the install target" ON)

set(SIMFOR_INC_DIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE STRING
        "Installation directory for include files, a relative path that will be
            joined with ${CMAKE_INSTALL_PREFIX} or an absolute path.")

file(READ include/simfor/internal/version.hpp version_h)
if (NOT version_h MATCHES "SIMFOR_VERSION ([0-9]+)([0-9][0-9])([0-9][0-9])")
    message(FATAL_ERROR "Cannot get SIMFOR_VERSION from version.h.")
endif ()
math(EXPR CPACK_PACKAGE_VERSION_MAJOR ${CMAKE_MATCH_1})
math(EXPR CPACK_PACKAGE_VERSION_MINOR ${CMAKE_MATCH_2})
math(EXPR CPACK_PACKAGE_VERSION_PATCH ${CMAKE_MATCH_3})
join(SIMFOR_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.
        ${CPACK_PACKAGE_VERSION_PATCH})
message(STATUS "{simfor} version: ${SIMFOR_VERSION}")

add_library(${PROJECT_NAME})

add_compile_options(-O3 -fPIC)

set(SIMFOR_INCLUDES
        # Public API includes
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/Gauss.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/GaussOmp.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/GaussMpi.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/Zeidel.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/ZeidelOmp.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/ZeidelMpi.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/LUdecomp.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/LUdecompOmp.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/SimpleIter.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/SimpleIterOmp.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/SimpleIterMpi.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/Gradients.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/GradientsOmp.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/GradientsMpi.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/Tridiagonal.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/TridiagonalOmp.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/TridiagonalMpi.hpp

        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/odu.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/elementary.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/classic.hpp

        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/plotter.h

        # Private API includes
        ${CMAKE_CURRENT_SOURCE_DIR}/include/simfor/internal/types.hpp
)

set(SIMFOR_SOURCES
        ${SIMFOR_INCLUDES}

        src/Gauss.cpp
        src/GaussOmp.cpp
        src/GaussMpi.cpp
        src/Zeidel.cpp
        src/ZeidelOmp.cpp
        src/ZeidelMpi.cpp
        src/LUdecomp.cpp
        src/LUdecompOmp.cpp
        src/SimpleIter.cpp
        src/SimpleIterOmp.cpp
        src/SimpleIterMpi.cpp
        src/Gradients.cpp
        src/GradientsOmp.cpp
        src/GradientsMpi.cpp
        src/Tridiagonal.cpp
        src/TridiagonalOmp.cpp
        src/TridiagonalMpi.cpp

        src/odu.cpp
        src/classic.cpp
        src/elementary.cpp

        src/plotter.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${SIMFOR_SOURCES})

target_sources(${PROJECT_NAME} PRIVATE ${SIMFOR_SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/${PROJECT_NAME}>
        $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
)

target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${SIMFOR_INC_DIR}>
)

target_include_directories(${PROJECT_NAME} PUBLIC
        ${Boost_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME} PUBLIC
        ${Boost_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        OpenMP::OpenMP_CXX
        GL
        glut
)

target_compile_options(${PROJECT_NAME} PUBLIC
        ${OpenMP_CXX_FLAGS}
)

if (SIMFOR_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif ()

if (SIMFOR_INSTALL)
    set(SIMFOR_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/simfor" CACHE STRING
            "Installation directory for cmake files, a relative path that will be
            joined with ${CMAKE_INSTALL_PREFIX} or an absolute path.")

    set(version_config ${PROJECT_BINARY_DIR}/simfor-config-version.cmake)
    set(project_config ${PROJECT_BINARY_DIR}/simfor-config.cmake)
    set(pkgconfig ${PROJECT_BINARY_DIR}/simfor.pc)
    set(targets_export_name simfor-targets)

    set(SIMFOR_LIB_DIR "${CMAKE_INSTALL_LIBDIR}" CACHE STRING
            "Installation directory for libraries, a relative path that will be joined
            to ${CMAKE_INSTALL_PREFIX} or an absolute path.")
    set(SIMFOR_PKGCONFIG_DIR "${CMAKE_INSTALL_LIBDIR}/pkgconfig" CACHE STRING
            "Installation directory for pkgconfig (.pc) files, a relative path that will
             be joined with ${CMAKE_INSTALL_PREFIX} or an absolute path.")

    set(SIMFOR_LIB_NAME ${PROJECT_NAME})

    write_basic_package_version_file(
            ${version_config}
            VERSION ${FMT_VERSION}
            COMPATIBILITY AnyNewerVersion
    )

    join_paths(libdir_for_pc_file "\${exec_prefix}" "${SIMFOR_LIB_DIR}")
    join_paths(includedir_for_pc_file "\${prefix}" "${SIMFOR_INC_DIR}")

    configure_file(
            "${PROJECT_SOURCE_DIR}/cmake/simfor.pc.in"
            "${pkgconfig}"
            @ONLY
    )
    configure_package_config_file(
            "${PROJECT_SOURCE_DIR}/cmake/simfor-config.cmake.in"
            "${project_config}"
            INSTALL_DESTINATION ${SIMFOR_CMAKE_DIR}
    )

    set(INSTALL_TARGETS simfor)

    # Install the library and headers.
    install(TARGETS ${INSTALL_TARGETS} EXPORT ${targets_export_name}
            LIBRARY DESTINATION ${SIMFOR_LIB_DIR}
            ARCHIVE DESTINATION ${SIMFOR_LIB_DIR}
            PUBLIC_HEADER DESTINATION "${SIMFOR_INC_DIR}/simfor"
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

    install(FILES ${SIMFOR_INCLUDES} DESTINATION include/simfor)

    # Use a namespace because CMake provides better diagnostics for namespaced
    # imported targets.
    export(TARGETS ${INSTALL_TARGETS} NAMESPACE simfor::
            FILE ${PROJECT_BINARY_DIR}/${targets_export_name}.cmake)

    # Install version, config and target files.
    install(
            FILES ${project_config} ${version_config}
            DESTINATION ${SIMFOR_CMAKE_DIR})
    install(EXPORT ${targets_export_name} DESTINATION ${SIMFOR_CMAKE_DIR}
            NAMESPACE simfor::)

    install(FILES "${pkgconfig}" DESTINATION "${SIMFOR_PKGCONFIG_DIR}")

endif ()

if (SIMFOR_BUILD_EXAMPLES)
    set(EXAMPLES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/examples")
    subdirlist(SUBDIRS ${EXAMPLES_DIR})
    foreach (subdir ${SUBDIRS})
        add_subdirectory("${EXAMPLES_DIR}/${subdir}")
    endforeach ()
endif ()

if (SIMFOR_BUILD_DOCS)
    include(doxygen)
endif ()