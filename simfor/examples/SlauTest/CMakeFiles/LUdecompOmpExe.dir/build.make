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
include examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/compiler_depend.make

# Include the progress variables for this target.
include examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/progress.make

# Include the compile flags for this target's objects.
include examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/flags.make

examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.o: examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/flags.make
examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.o: ../examples/SlauTest/LUdecompOmpExe.cpp
examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.o: examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vadim/programs/cpp/simfor/simfor/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.o"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/SlauTest && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.o -MF CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.o.d -o CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.o -c /home/vadim/programs/cpp/simfor/examples/SlauTest/LUdecompOmpExe.cpp

examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.i"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/SlauTest && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vadim/programs/cpp/simfor/examples/SlauTest/LUdecompOmpExe.cpp > CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.i

examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.s"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/SlauTest && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vadim/programs/cpp/simfor/examples/SlauTest/LUdecompOmpExe.cpp -o CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.s

# Object files for target LUdecompOmpExe
LUdecompOmpExe_OBJECTS = \
"CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.o"

# External object files for target LUdecompOmpExe
LUdecompOmpExe_EXTERNAL_OBJECTS =

examples/SlauTest/LUdecompOmpExe: examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/LUdecompOmpExe.cpp.o
examples/SlauTest/LUdecompOmpExe: examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/build.make
examples/SlauTest/LUdecompOmpExe: libsimfor.a
examples/SlauTest/LUdecompOmpExe: /usr/local/lib/libboost_mpi.so.1.84.0
examples/SlauTest/LUdecompOmpExe: /usr/local/lib/libboost_serialization.so.1.84.0
examples/SlauTest/LUdecompOmpExe: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
examples/SlauTest/LUdecompOmpExe: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
examples/SlauTest/LUdecompOmpExe: /usr/lib/x86_64-linux-gnu/libOpenGL.so
examples/SlauTest/LUdecompOmpExe: /usr/lib/x86_64-linux-gnu/libGLX.so
examples/SlauTest/LUdecompOmpExe: /usr/lib/x86_64-linux-gnu/libGLU.so
examples/SlauTest/LUdecompOmpExe: /usr/lib/x86_64-linux-gnu/libglut.so
examples/SlauTest/LUdecompOmpExe: /usr/lib/gcc/x86_64-linux-gnu/11/libgomp.so
examples/SlauTest/LUdecompOmpExe: /usr/lib/x86_64-linux-gnu/libpthread.a
examples/SlauTest/LUdecompOmpExe: examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/vadim/programs/cpp/simfor/simfor/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable LUdecompOmpExe"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/SlauTest && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LUdecompOmpExe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/build: examples/SlauTest/LUdecompOmpExe
.PHONY : examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/build

examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/clean:
	cd /home/vadim/programs/cpp/simfor/simfor/examples/SlauTest && $(CMAKE_COMMAND) -P CMakeFiles/LUdecompOmpExe.dir/cmake_clean.cmake
.PHONY : examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/clean

examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/depend:
	cd /home/vadim/programs/cpp/simfor/simfor && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vadim/programs/cpp/simfor /home/vadim/programs/cpp/simfor/examples/SlauTest /home/vadim/programs/cpp/simfor/simfor /home/vadim/programs/cpp/simfor/simfor/examples/SlauTest /home/vadim/programs/cpp/simfor/simfor/examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/SlauTest/CMakeFiles/LUdecompOmpExe.dir/depend

