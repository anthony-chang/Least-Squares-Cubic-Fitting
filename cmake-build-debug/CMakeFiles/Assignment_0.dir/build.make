# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.10

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2017.3\bin\cmake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2017.3\bin\cmake\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "C:\Users\Anthony\Desktop\School\Grade 12\Semester 1\AP Physics C-Mr. van Bemmel\Assignment 0"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "C:\Users\Anthony\Desktop\School\Grade 12\Semester 1\AP Physics C-Mr. van Bemmel\Assignment 0\cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/Assignment_0.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Assignment_0.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Assignment_0.dir/flags.make

CMakeFiles/Assignment_0.dir/main.cpp.obj: CMakeFiles/Assignment_0.dir/flags.make
CMakeFiles/Assignment_0.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\Anthony\Desktop\School\Grade 12\Semester 1\AP Physics C-Mr. van Bemmel\Assignment 0\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Assignment_0.dir/main.cpp.obj"
	C:\PROGRA~2\MINGW-~1\I686-7~1.0-P\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\Assignment_0.dir\main.cpp.obj -c "C:\Users\Anthony\Desktop\School\Grade 12\Semester 1\AP Physics C-Mr. van Bemmel\Assignment 0\main.cpp"

CMakeFiles/Assignment_0.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Assignment_0.dir/main.cpp.i"
	C:\PROGRA~2\MINGW-~1\I686-7~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Users\Anthony\Desktop\School\Grade 12\Semester 1\AP Physics C-Mr. van Bemmel\Assignment 0\main.cpp" > CMakeFiles\Assignment_0.dir\main.cpp.i

CMakeFiles/Assignment_0.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Assignment_0.dir/main.cpp.s"
	C:\PROGRA~2\MINGW-~1\I686-7~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Users\Anthony\Desktop\School\Grade 12\Semester 1\AP Physics C-Mr. van Bemmel\Assignment 0\main.cpp" -o CMakeFiles\Assignment_0.dir\main.cpp.s

CMakeFiles/Assignment_0.dir/main.cpp.obj.requires:

.PHONY : CMakeFiles/Assignment_0.dir/main.cpp.obj.requires

CMakeFiles/Assignment_0.dir/main.cpp.obj.provides: CMakeFiles/Assignment_0.dir/main.cpp.obj.requires
	$(MAKE) -f CMakeFiles\Assignment_0.dir\build.make CMakeFiles/Assignment_0.dir/main.cpp.obj.provides.build
.PHONY : CMakeFiles/Assignment_0.dir/main.cpp.obj.provides

CMakeFiles/Assignment_0.dir/main.cpp.obj.provides.build: CMakeFiles/Assignment_0.dir/main.cpp.obj


# Object files for target Assignment_0
Assignment_0_OBJECTS = \
"CMakeFiles/Assignment_0.dir/main.cpp.obj"

# External object files for target Assignment_0
Assignment_0_EXTERNAL_OBJECTS =

Assignment_0.exe: CMakeFiles/Assignment_0.dir/main.cpp.obj
Assignment_0.exe: CMakeFiles/Assignment_0.dir/build.make
Assignment_0.exe: CMakeFiles/Assignment_0.dir/linklibs.rsp
Assignment_0.exe: CMakeFiles/Assignment_0.dir/objects1.rsp
Assignment_0.exe: CMakeFiles/Assignment_0.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\Users\Anthony\Desktop\School\Grade 12\Semester 1\AP Physics C-Mr. van Bemmel\Assignment 0\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Assignment_0.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\Assignment_0.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Assignment_0.dir/build: Assignment_0.exe

.PHONY : CMakeFiles/Assignment_0.dir/build

CMakeFiles/Assignment_0.dir/requires: CMakeFiles/Assignment_0.dir/main.cpp.obj.requires

.PHONY : CMakeFiles/Assignment_0.dir/requires

CMakeFiles/Assignment_0.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\Assignment_0.dir\cmake_clean.cmake
.PHONY : CMakeFiles/Assignment_0.dir/clean

CMakeFiles/Assignment_0.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\Users\Anthony\Desktop\School\Grade 12\Semester 1\AP Physics C-Mr. van Bemmel\Assignment 0" "C:\Users\Anthony\Desktop\School\Grade 12\Semester 1\AP Physics C-Mr. van Bemmel\Assignment 0" "C:\Users\Anthony\Desktop\School\Grade 12\Semester 1\AP Physics C-Mr. van Bemmel\Assignment 0\cmake-build-debug" "C:\Users\Anthony\Desktop\School\Grade 12\Semester 1\AP Physics C-Mr. van Bemmel\Assignment 0\cmake-build-debug" "C:\Users\Anthony\Desktop\School\Grade 12\Semester 1\AP Physics C-Mr. van Bemmel\Assignment 0\cmake-build-debug\CMakeFiles\Assignment_0.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/Assignment_0.dir/depend
