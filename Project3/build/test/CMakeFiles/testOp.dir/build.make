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
include test/CMakeFiles/testOp.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/testOp.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/testOp.dir/flags.make

test/CMakeFiles/testOp.dir/testOperator.cpp.o: test/CMakeFiles/testOp.dir/flags.make
test/CMakeFiles/testOp.dir/testOperator.cpp.o: ../test/testOperator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/shuang/numPDEs/Project3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/testOp.dir/testOperator.cpp.o"
	cd /home/shuang/numPDEs/Project3/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/testOp.dir/testOperator.cpp.o -c /home/shuang/numPDEs/Project3/test/testOperator.cpp

test/CMakeFiles/testOp.dir/testOperator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testOp.dir/testOperator.cpp.i"
	cd /home/shuang/numPDEs/Project3/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/shuang/numPDEs/Project3/test/testOperator.cpp > CMakeFiles/testOp.dir/testOperator.cpp.i

test/CMakeFiles/testOp.dir/testOperator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testOp.dir/testOperator.cpp.s"
	cd /home/shuang/numPDEs/Project3/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/shuang/numPDEs/Project3/test/testOperator.cpp -o CMakeFiles/testOp.dir/testOperator.cpp.s

# Object files for target testOp
testOp_OBJECTS = \
"CMakeFiles/testOp.dir/testOperator.cpp.o"

# External object files for target testOp
testOp_EXTERNAL_OBJECTS =

test/testOp: test/CMakeFiles/testOp.dir/testOperator.cpp.o
test/testOp: test/CMakeFiles/testOp.dir/build.make
test/testOp: src/libMultiGrid.a
test/testOp: test/CMakeFiles/testOp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/shuang/numPDEs/Project3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable testOp"
	cd /home/shuang/numPDEs/Project3/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testOp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/testOp.dir/build: test/testOp

.PHONY : test/CMakeFiles/testOp.dir/build

test/CMakeFiles/testOp.dir/clean:
	cd /home/shuang/numPDEs/Project3/build/test && $(CMAKE_COMMAND) -P CMakeFiles/testOp.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/testOp.dir/clean

test/CMakeFiles/testOp.dir/depend:
	cd /home/shuang/numPDEs/Project3/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shuang/numPDEs/Project3 /home/shuang/numPDEs/Project3/test /home/shuang/numPDEs/Project3/build /home/shuang/numPDEs/Project3/build/test /home/shuang/numPDEs/Project3/build/test/CMakeFiles/testOp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/testOp.dir/depend

