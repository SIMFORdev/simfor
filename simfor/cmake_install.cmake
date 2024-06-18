# Install script for directory: /home/vadim/programs/cpp/simfor

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/vadim/programs/cpp/simfor/simfor/libsimfor.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simfor" TYPE FILE FILES
    "/home/vadim/programs/cpp/simfor/include/simfor/Gauss.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/GaussOmp.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/GaussMpi.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/Zeidel.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/ZeidelOmp.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/ZeidelMpi.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/LUdecomp.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/LUdecompOmp.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/SimpleIter.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/SimpleIterOmp.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/SimpleIterMpi.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/Gradients.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/GradientsOmp.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/GradientsMpi.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/Tridiagonal.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/TridiagonalOmp.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/TridiagonalMpi.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/odu.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/elementary.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/classic.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/classic_omp.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/strassen.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/strassen_omp.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/strassen_recursive.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/strassen_recursive_omp.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/cblock.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/cblock_omp.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/multiplyVectorMatrix.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/multMatrVec.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/plotter.h"
    "/home/vadim/programs/cpp/simfor/include/simfor/types.hpp"
    "/home/vadim/programs/cpp/simfor/include/simfor/version.hpp"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/simfor" TYPE FILE FILES
    "/home/vadim/programs/cpp/simfor/simfor/simfor-config.cmake"
    "/home/vadim/programs/cpp/simfor/simfor/simfor-config-version.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/simfor/simfor-targets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/simfor/simfor-targets.cmake"
         "/home/vadim/programs/cpp/simfor/simfor/CMakeFiles/Export/lib/cmake/simfor/simfor-targets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/simfor/simfor-targets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/simfor/simfor-targets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/simfor" TYPE FILE FILES "/home/vadim/programs/cpp/simfor/simfor/CMakeFiles/Export/lib/cmake/simfor/simfor-targets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/simfor" TYPE FILE FILES "/home/vadim/programs/cpp/simfor/simfor/CMakeFiles/Export/lib/cmake/simfor/simfor-targets-noconfig.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/vadim/programs/cpp/simfor/simfor/simfor.pc")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/vadim/programs/cpp/simfor/simfor/examples/BOOST_ODE/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/vadim/programs/cpp/simfor/simfor/examples/MatrixTest/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/vadim/programs/cpp/simfor/simfor/examples/ODE/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/vadim/programs/cpp/simfor/simfor/examples/SlauTest/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/vadim/programs/cpp/simfor/simfor/examples/plotter/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/vadim/programs/cpp/simfor/simfor/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
