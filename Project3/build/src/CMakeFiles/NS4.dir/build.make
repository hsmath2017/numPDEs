# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/shuang/numPDEs/Project3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/shuang/numPDEs/Project3/build

# Include any dependencies generated for this target.
include src/CMakeFiles/NS4.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/NS4.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/NS4.dir/flags.make

src/CMakeFiles/NS4.dir/Core/Config.cpp.o: src/CMakeFiles/NS4.dir/flags.make
src/CMakeFiles/NS4.dir/Core/Config.cpp.o: ../src/Core/Config.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/shuang/numPDEs/Project3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/NS4.dir/Core/Config.cpp.o"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NS4.dir/Core/Config.cpp.o -c /home/shuang/numPDEs/Project3/src/Core/Config.cpp

src/CMakeFiles/NS4.dir/Core/Config.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NS4.dir/Core/Config.cpp.i"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/shuang/numPDEs/Project3/src/Core/Config.cpp > CMakeFiles/NS4.dir/Core/Config.cpp.i

src/CMakeFiles/NS4.dir/Core/Config.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NS4.dir/Core/Config.cpp.s"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/shuang/numPDEs/Project3/src/Core/Config.cpp -o CMakeFiles/NS4.dir/Core/Config.cpp.s

src/CMakeFiles/NS4.dir/Core/numlib.cpp.o: src/CMakeFiles/NS4.dir/flags.make
src/CMakeFiles/NS4.dir/Core/numlib.cpp.o: ../src/Core/numlib.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/shuang/numPDEs/Project3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/NS4.dir/Core/numlib.cpp.o"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NS4.dir/Core/numlib.cpp.o -c /home/shuang/numPDEs/Project3/src/Core/numlib.cpp

src/CMakeFiles/NS4.dir/Core/numlib.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NS4.dir/Core/numlib.cpp.i"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/shuang/numPDEs/Project3/src/Core/numlib.cpp > CMakeFiles/NS4.dir/Core/numlib.cpp.i

src/CMakeFiles/NS4.dir/Core/numlib.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NS4.dir/Core/numlib.cpp.s"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/shuang/numPDEs/Project3/src/Core/numlib.cpp -o CMakeFiles/NS4.dir/Core/numlib.cpp.s

src/CMakeFiles/NS4.dir/Core/Curve.cpp.o: src/CMakeFiles/NS4.dir/flags.make
src/CMakeFiles/NS4.dir/Core/Curve.cpp.o: ../src/Core/Curve.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/shuang/numPDEs/Project3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/NS4.dir/Core/Curve.cpp.o"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NS4.dir/Core/Curve.cpp.o -c /home/shuang/numPDEs/Project3/src/Core/Curve.cpp

src/CMakeFiles/NS4.dir/Core/Curve.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NS4.dir/Core/Curve.cpp.i"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/shuang/numPDEs/Project3/src/Core/Curve.cpp > CMakeFiles/NS4.dir/Core/Curve.cpp.i

src/CMakeFiles/NS4.dir/Core/Curve.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NS4.dir/Core/Curve.cpp.s"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/shuang/numPDEs/Project3/src/Core/Curve.cpp -o CMakeFiles/NS4.dir/Core/Curve.cpp.s

src/CMakeFiles/NS4.dir/Core/Wrapper_OpenMP.cpp.o: src/CMakeFiles/NS4.dir/flags.make
src/CMakeFiles/NS4.dir/Core/Wrapper_OpenMP.cpp.o: ../src/Core/Wrapper_OpenMP.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/shuang/numPDEs/Project3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/CMakeFiles/NS4.dir/Core/Wrapper_OpenMP.cpp.o"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NS4.dir/Core/Wrapper_OpenMP.cpp.o -c /home/shuang/numPDEs/Project3/src/Core/Wrapper_OpenMP.cpp

src/CMakeFiles/NS4.dir/Core/Wrapper_OpenMP.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NS4.dir/Core/Wrapper_OpenMP.cpp.i"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/shuang/numPDEs/Project3/src/Core/Wrapper_OpenMP.cpp > CMakeFiles/NS4.dir/Core/Wrapper_OpenMP.cpp.i

src/CMakeFiles/NS4.dir/Core/Wrapper_OpenMP.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NS4.dir/Core/Wrapper_OpenMP.cpp.s"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/shuang/numPDEs/Project3/src/Core/Wrapper_OpenMP.cpp -o CMakeFiles/NS4.dir/Core/Wrapper_OpenMP.cpp.s

src/CMakeFiles/NS4.dir/MultiGrid/MGSolver.cpp.o: src/CMakeFiles/NS4.dir/flags.make
src/CMakeFiles/NS4.dir/MultiGrid/MGSolver.cpp.o: ../src/MultiGrid/MGSolver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/shuang/numPDEs/Project3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/CMakeFiles/NS4.dir/MultiGrid/MGSolver.cpp.o"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NS4.dir/MultiGrid/MGSolver.cpp.o -c /home/shuang/numPDEs/Project3/src/MultiGrid/MGSolver.cpp

src/CMakeFiles/NS4.dir/MultiGrid/MGSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NS4.dir/MultiGrid/MGSolver.cpp.i"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/shuang/numPDEs/Project3/src/MultiGrid/MGSolver.cpp > CMakeFiles/NS4.dir/MultiGrid/MGSolver.cpp.i

src/CMakeFiles/NS4.dir/MultiGrid/MGSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NS4.dir/MultiGrid/MGSolver.cpp.s"
	cd /home/shuang/numPDEs/Project3/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/shuang/numPDEs/Project3/src/MultiGrid/MGSolver.cpp -o CMakeFiles/NS4.dir/MultiGrid/MGSolver.cpp.s

# Object files for target NS4
NS4_OBJECTS = \
"CMakeFiles/NS4.dir/Core/Config.cpp.o" \
"CMakeFiles/NS4.dir/Core/numlib.cpp.o" \
"CMakeFiles/NS4.dir/Core/Curve.cpp.o" \
"CMakeFiles/NS4.dir/Core/Wrapper_OpenMP.cpp.o" \
"CMakeFiles/NS4.dir/MultiGrid/MGSolver.cpp.o"

# External object files for target NS4
NS4_EXTERNAL_OBJECTS =

src/libNS4.a: src/CMakeFiles/NS4.dir/Core/Config.cpp.o
src/libNS4.a: src/CMakeFiles/NS4.dir/Core/numlib.cpp.o
src/libNS4.a: src/CMakeFiles/NS4.dir/Core/Curve.cpp.o
src/libNS4.a: src/CMakeFiles/NS4.dir/Core/Wrapper_OpenMP.cpp.o
src/libNS4.a: src/CMakeFiles/NS4.dir/MultiGrid/MGSolver.cpp.o
src/libNS4.a: src/CMakeFiles/NS4.dir/build.make
src/libNS4.a: src/CMakeFiles/NS4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/shuang/numPDEs/Project3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX static library libNS4.a"
	cd /home/shuang/numPDEs/Project3/build/src && $(CMAKE_COMMAND) -P CMakeFiles/NS4.dir/cmake_clean_target.cmake
	cd /home/shuang/numPDEs/Project3/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/NS4.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/NS4.dir/build: src/libNS4.a

.PHONY : src/CMakeFiles/NS4.dir/build

src/CMakeFiles/NS4.dir/clean:
	cd /home/shuang/numPDEs/Project3/build/src && $(CMAKE_COMMAND) -P CMakeFiles/NS4.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/NS4.dir/clean

src/CMakeFiles/NS4.dir/depend:
	cd /home/shuang/numPDEs/Project3/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shuang/numPDEs/Project3 /home/shuang/numPDEs/Project3/src /home/shuang/numPDEs/Project3/build /home/shuang/numPDEs/Project3/build/src /home/shuang/numPDEs/Project3/build/src/CMakeFiles/NS4.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/NS4.dir/depend

