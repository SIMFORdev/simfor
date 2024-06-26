find_package(Doxygen REQUIRED)
if (DOXYGEN_FOUND)
    set(SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/src
            ${CMAKE_CURRENT_SOURCE_DIR}/include)
    set(DOXYGEN_PROJECT_NAME "simfor")
    set(DOXYGEN_PROJECT_NUMBER "0.0.2")
    set(DOXYGEN_GENERATE_MAN "YES")
    set(DOXYGEN_OPTIMIZE_OUTPUT_FOR_C "YES")
    set(DOXYGEN_EXTRACT_ALL "YES")
    set(DOXYGEN_EXAMPLE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/examples)
    set(DOXYGEN_EXAMPLE_RECURSIVE "YES")
    set(DOXYGEN_GENERATE_TREEVIEW "YES")

    doxygen_add_docs(txt ${SOURCES})

    set(DOXYGEN_INPUT ${PROJECT_BINARY_DIR}/Doxyfile.txt)
    set(DOXYGEN_OUTPUT /doc/html)

    add_custom_command(
            OUTPUT ${DOXYGEN_OUTPUT}
            COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
            COMMAND ${CMAKE_COMMAND} -E echo "Done."
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            DEPENDS ${DOXYGEN_INPUT}
    )

    add_custom_target(apidoc ALL DEPENDS ${DOXYGEN_OUTPUT})

endif (DOXYGEN_FOUND)
