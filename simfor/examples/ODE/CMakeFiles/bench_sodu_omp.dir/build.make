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
include examples/ODE/CMakeFiles/bench_sodu_omp.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include examples/ODE/CMakeFiles/bench_sodu_omp.dir/compiler_depend.make

# Include the progress variables for this target.
include examples/ODE/CMakeFiles/bench_sodu_omp.dir/progress.make

# Include the compile flags for this target's objects.
include examples/ODE/CMakeFiles/bench_sodu_omp.dir/flags.make

examples/ODE/CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.o: examples/ODE/CMakeFiles/bench_sodu_omp.dir/flags.make
examples/ODE/CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.o: ../examples/ODE/bench_sodu_omp.cpp
examples/ODE/CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.o: examples/ODE/CMakeFiles/bench_sodu_omp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vadim/programs/cpp/simfor/simfor/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/ODE/CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.o"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/ODE && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT examples/ODE/CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.o -MF CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.o.d -o CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.o -c /home/vadim/programs/cpp/simfor/examples/ODE/bench_sodu_omp.cpp

examples/ODE/CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.i"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/ODE && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vadim/programs/cpp/simfor/examples/ODE/bench_sodu_omp.cpp > CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.i

examples/ODE/CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.s"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/ODE && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vadim/programs/cpp/simfor/examples/ODE/bench_sodu_omp.cpp -o CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.s

examples/ODE/CMakeFiles/bench_sodu_omp.dir/sodu.cpp.o: examples/ODE/CMakeFiles/bench_sodu_omp.dir/flags.make
examples/ODE/CMakeFiles/bench_sodu_omp.dir/sodu.cpp.o: ../examples/ODE/sodu.cpp
examples/ODE/CMakeFiles/bench_sodu_omp.dir/sodu.cpp.o: examples/ODE/CMakeFiles/bench_sodu_omp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vadim/programs/cpp/simfor/simfor/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object examples/ODE/CMakeFiles/bench_sodu_omp.dir/sodu.cpp.o"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/ODE && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT examples/ODE/CMakeFiles/bench_sodu_omp.dir/sodu.cpp.o -MF CMakeFiles/bench_sodu_omp.dir/sodu.cpp.o.d -o CMakeFiles/bench_sodu_omp.dir/sodu.cpp.o -c /home/vadim/programs/cpp/simfor/examples/ODE/sodu.cpp

examples/ODE/CMakeFiles/bench_sodu_omp.dir/sodu.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bench_sodu_omp.dir/sodu.cpp.i"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/ODE && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vadim/programs/cpp/simfor/examples/ODE/sodu.cpp > CMakeFiles/bench_sodu_omp.dir/sodu.cpp.i

examples/ODE/CMakeFiles/bench_sodu_omp.dir/sodu.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bench_sodu_omp.dir/sodu.cpp.s"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/ODE && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vadim/programs/cpp/simfor/examples/ODE/sodu.cpp -o CMakeFiles/bench_sodu_omp.dir/sodu.cpp.s

# Object files for target bench_sodu_omp
bench_sodu_omp_OBJECTS = \
"CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.o" \
"CMakeFiles/bench_sodu_omp.dir/sodu.cpp.o"

# External object files for target bench_sodu_omp
bench_sodu_omp_EXTERNAL_OBJECTS =

examples/ODE/bench_sodu_omp: examples/ODE/CMakeFiles/bench_sodu_omp.dir/bench_sodu_omp.cpp.o
examples/ODE/bench_sodu_omp: examples/ODE/CMakeFiles/bench_sodu_omp.dir/sodu.cpp.o
examples/ODE/bench_sodu_omp: examples/ODE/CMakeFiles/bench_sodu_omp.dir/build.make
examples/ODE/bench_sodu_omp: libsimfor.a
examples/ODE/bench_sodu_omp: /usr/local/lib/libboost_mpi.so.1.84.0
examples/ODE/bench_sodu_omp: /usr/local/lib/libboost_serialization.so.1.84.0
examples/ODE/bench_sodu_omp: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
examples/ODE/bench_sodu_omp: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
examples/ODE/bench_sodu_omp: /usr/lib/x86_64-linux-gnu/libOpenGL.so
examples/ODE/bench_sodu_omp: /usr/lib/x86_64-linux-gnu/libGLX.so
examples/ODE/bench_sodu_omp: /usr/lib/x86_64-linux-gnu/libGLU.so
examples/ODE/bench_sodu_omp: /usr/lib/x86_64-linux-gnu/libglut.so
examples/ODE/bench_sodu_omp: /usr/lib/gcc/x86_64-linux-gnu/11/libgomp.so
examples/ODE/bench_sodu_omp: /usr/lib/x86_64-linux-gnu/libpthread.a
examples/ODE/bench_sodu_omp: examples/ODE/CMakeFiles/bench_sodu_omp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/vadim/programs/cpp/simfor/simfor/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable bench_sodu_omp"
	cd /home/vadim/programs/cpp/simfor/simfor/examples/ODE && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bench_sodu_omp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/ODE/CMakeFiles/bench_sodu_omp.dir/build: examples/ODE/bench_sodu_omp
.PHONY : examples/ODE/CMakeFiles/bench_sodu_omp.dir/build

examples/ODE/CMakeFiles/bench_sodu_omp.dir/clean:
	cd /home/vadim/programs/cpp/simfor/simfor/examples/ODE && $(CMAKE_COMMAND) -P CMakeFiles/bench_sodu_omp.dir/cmake_clean.cmake
.PHONY : examples/ODE/CMakeFiles/bench_sodu_omp.dir/clean

examples/ODE/CMakeFiles/bench_sodu_omp.dir/depend:
	cd /home/vadim/programs/cpp/simfor/simfor && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vadim/programs/cpp/simfor /home/vadim/programs/cpp/simfor/examples/ODE /home/vadim/programs/cpp/simfor/simfor /home/vadim/programs/cpp/simfor/simfor/examples/ODE /home/vadim/programs/cpp/simfor/simfor/examples/ODE/CMakeFiles/bench_sodu_omp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/ODE/CMakeFiles/bench_sodu_omp.dir/depend

