# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/vadim/programs/cpp/simfor

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/vadim/programs/cpp/simfor/simfor

# Include any dependencies generated for this target.
include examples/plotter/CMakeFiles/test_2d.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include examples/plotter/CMakeFiles/test_2d.dir/compiler_depend.make

# Include the progress variables for this target.
include examples/plotter/CMakeFiles/test_2d.dir/progress.make

# Include the compile flags for this target's objects.
include examples/plotter/CMakeFiles/test_2d.dir/flags.make

examples/plotter/CMakeFiles/test_2d.dir/test_2d.cpp.o: examples/plotter/CMakeFiles/test_2d.dir/flags.make
examples/plotter/CMakeFiles/test_2d.dir/test_2d.cpp.o: ../examples/plotter/test_2d.cpp
examples/plotter/CMakeFiles/test_2d.dir/test_2d.cpp.o: examples/plotter/CMakeFiles/test_2d.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vadim/programs/cpp/simfor/simfor/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/plotter/CMakeFiles/test_2d.dir/test_2d.cpp.o"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/plotter && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT examples/plotter/CMakeFiles/test_2d.dir/test_2d.cpp.o -MF CMakeFiles/test_2d.dir/test_2d.cpp.o.d -o CMakeFiles/test_2d.dir/test_2d.cpp.o -c /home/vadim/programs/cpp/simfor/examples/plotter/test_2d.cpp

examples/plotter/CMakeFiles/test_2d.dir/test_2d.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_2d.dir/test_2d.cpp.i"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/plotter && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vadim/programs/cpp/simfor/examples/plotter/test_2d.cpp > CMakeFiles/test_2d.dir/test_2d.cpp.i

examples/plotter/CMakeFiles/test_2d.dir/test_2d.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_2d.dir/test_2d.cpp.s"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/plotter && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vadim/programs/cpp/simfor/examples/plotter/test_2d.cpp -o CMakeFiles/test_2d.dir/test_2d.cpp.s

# Object files for target test_2d
test_2d_OBJECTS = \
"CMakeFiles/test_2d.dir/test_2d.cpp.o"

# External object files for target test_2d
test_2d_EXTERNAL_OBJECTS =

examples/plotter/test_2d: examples/plotter/CMakeFiles/test_2d.dir/test_2d.cpp.o
examples/plotter/test_2d: examples/plotter/CMakeFiles/test_2d.dir/build.make
examples/plotter/test_2d: libsimfor.a
examples/plotter/test_2d: /usr/local/lib/libboost_mpi.so.1.84.0
examples/plotter/test_2d: /usr/local/lib/libboost_serialization.so.1.84.0
examples/plotter/test_2d: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
examples/plotter/test_2d: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
examples/plotter/test_2d: /usr/lib/x86_64-linux-gnu/libOpenGL.so
examples/plotter/test_2d: /usr/lib/x86_64-linux-gnu/libGLX.so
examples/plotter/test_2d: /usr/lib/x86_64-linux-gnu/libGLU.so
examples/plotter/test_2d: /usr/lib/x86_64-linux-gnu/libglut.so
examples/plotter/test_2d: /usr/lib/gcc/x86_64-linux-gnu/11/libgomp.so
examples/plotter/test_2d: /usr/lib/x86_64-linux-gnu/libpthread.a
examples/plotter/test_2d: examples/plotter/CMakeFiles/test_2d.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/vadim/programs/cpp/simfor/simfor/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_2d"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/plotter && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_2d.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/plotter/CMakeFiles/test_2d.dir/build: examples/plotter/test_2d
.PHONY : examples/plotter/CMakeFiles/test_2d.dir/build

examples/plotter/CMakeFiles/test_2d.dir/clean:
	cd /home/vadim/programs/cpp/simfor/simfor/examples/plotter && $(CMAKE_COMMAND) -P CMakeFiles/test_2d.dir/cmake_clean.cmake
.PHONY : examples/plotter/CMakeFiles/test_2d.dir/clean

examples/plotter/CMakeFiles/test_2d.dir/depend:
	cd /home/vadim/programs/cpp/simfor/simfor && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vadim/programs/cpp/simfor /home/vadim/programs/cpp/simfor/examples/plotter /home/vadim/programs/cpp/simfor/simfor /home/vadim/programs/cpp/simfor/simfor/examples/plotter /home/vadim/programs/cpp/simfor/simfor/examples/plotter/CMakeFiles/test_2d.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/plotter/CMakeFiles/test_2d.dir/depend

