# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

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
CMAKE_SOURCE_DIR = /home/andrej/studyspace/courseProject/pclTest

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/andrej/studyspace/courseProject/pclTest

# Include any dependencies generated for this target.
include CMakeFiles/openniViewer.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/openniViewer.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/openniViewer.dir/flags.make

CMakeFiles/openniViewer.dir/src/main.cpp.o: CMakeFiles/openniViewer.dir/flags.make
CMakeFiles/openniViewer.dir/src/main.cpp.o: src/main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/andrej/studyspace/courseProject/pclTest/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/openniViewer.dir/src/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/openniViewer.dir/src/main.cpp.o -c /home/andrej/studyspace/courseProject/pclTest/src/main.cpp

CMakeFiles/openniViewer.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/openniViewer.dir/src/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/andrej/studyspace/courseProject/pclTest/src/main.cpp > CMakeFiles/openniViewer.dir/src/main.cpp.i

CMakeFiles/openniViewer.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/openniViewer.dir/src/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/andrej/studyspace/courseProject/pclTest/src/main.cpp -o CMakeFiles/openniViewer.dir/src/main.cpp.s

CMakeFiles/openniViewer.dir/src/main.cpp.o.requires:
.PHONY : CMakeFiles/openniViewer.dir/src/main.cpp.o.requires

CMakeFiles/openniViewer.dir/src/main.cpp.o.provides: CMakeFiles/openniViewer.dir/src/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/openniViewer.dir/build.make CMakeFiles/openniViewer.dir/src/main.cpp.o.provides.build
.PHONY : CMakeFiles/openniViewer.dir/src/main.cpp.o.provides

CMakeFiles/openniViewer.dir/src/main.cpp.o.provides.build: CMakeFiles/openniViewer.dir/src/main.cpp.o

# Object files for target openniViewer
openniViewer_OBJECTS = \
"CMakeFiles/openniViewer.dir/src/main.cpp.o"

# External object files for target openniViewer
openniViewer_EXTERNAL_OBJECTS =

openniViewer: CMakeFiles/openniViewer.dir/src/main.cpp.o
openniViewer: CMakeFiles/openniViewer.dir/build.make
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_system.so
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_thread.so
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_serialization.so
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_chrono.so
openniViewer: /usr/lib/x86_64-linux-gnu/libpthread.so
openniViewer: /usr/local/lib/libpcl_common.so
openniViewer: /usr/local/lib/libpcl_octree.so
openniViewer: /usr/lib/libOpenNI.so
openniViewer: /usr/lib/libOpenNI2.so
openniViewer: /usr/local/lib/libpcl_io.so
openniViewer: /usr/lib/x86_64-linux-gnu/libflann_cpp_s.a
openniViewer: /usr/local/lib/libpcl_kdtree.so
openniViewer: /usr/local/lib/libpcl_search.so
openniViewer: /usr/local/lib/libpcl_sample_consensus.so
openniViewer: /usr/local/lib/libpcl_filters.so
openniViewer: /usr/local/lib/libpcl_features.so
openniViewer: /usr/local/lib/libpcl_ml.so
openniViewer: /usr/local/lib/libpcl_segmentation.so
openniViewer: /usr/lib/x86_64-linux-gnu/libqhull.so
openniViewer: /usr/local/lib/libpcl_surface.so
openniViewer: /usr/local/lib/libpcl_registration.so
openniViewer: /usr/local/lib/libpcl_recognition.so
openniViewer: /usr/local/lib/libpcl_keypoints.so
openniViewer: /usr/local/lib/libpcl_visualization.so
openniViewer: /usr/local/lib/libpcl_tracking.so
openniViewer: /usr/local/lib/libpcl_stereo.so
openniViewer: /usr/local/lib/libpcl_outofcore.so
openniViewer: /usr/local/lib/libpcl_people.so
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_system.so
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_thread.so
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_serialization.so
openniViewer: /usr/lib/x86_64-linux-gnu/libboost_chrono.so
openniViewer: /usr/lib/x86_64-linux-gnu/libpthread.so
openniViewer: /usr/lib/x86_64-linux-gnu/libqhull.so
openniViewer: /usr/lib/libOpenNI.so
openniViewer: /usr/lib/libOpenNI2.so
openniViewer: /usr/lib/x86_64-linux-gnu/libflann_cpp_s.a
openniViewer: /usr/lib/libvtkGenericFiltering.so.5.8.0
openniViewer: /usr/lib/libvtkGeovis.so.5.8.0
openniViewer: /usr/lib/libvtkCharts.so.5.8.0
openniViewer: /usr/local/lib/libpcl_common.so
openniViewer: /usr/local/lib/libpcl_octree.so
openniViewer: /usr/local/lib/libpcl_io.so
openniViewer: /usr/local/lib/libpcl_kdtree.so
openniViewer: /usr/local/lib/libpcl_search.so
openniViewer: /usr/local/lib/libpcl_sample_consensus.so
openniViewer: /usr/local/lib/libpcl_filters.so
openniViewer: /usr/local/lib/libpcl_features.so
openniViewer: /usr/local/lib/libpcl_ml.so
openniViewer: /usr/local/lib/libpcl_segmentation.so
openniViewer: /usr/local/lib/libpcl_surface.so
openniViewer: /usr/local/lib/libpcl_registration.so
openniViewer: /usr/local/lib/libpcl_recognition.so
openniViewer: /usr/local/lib/libpcl_keypoints.so
openniViewer: /usr/local/lib/libpcl_visualization.so
openniViewer: /usr/local/lib/libpcl_tracking.so
openniViewer: /usr/local/lib/libpcl_stereo.so
openniViewer: /usr/local/lib/libpcl_outofcore.so
openniViewer: /usr/local/lib/libpcl_people.so
openniViewer: /usr/lib/libvtkViews.so.5.8.0
openniViewer: /usr/lib/libvtkInfovis.so.5.8.0
openniViewer: /usr/lib/libvtkWidgets.so.5.8.0
openniViewer: /usr/lib/libvtkVolumeRendering.so.5.8.0
openniViewer: /usr/lib/libvtkHybrid.so.5.8.0
openniViewer: /usr/lib/libvtkParallel.so.5.8.0
openniViewer: /usr/lib/libvtkRendering.so.5.8.0
openniViewer: /usr/lib/libvtkImaging.so.5.8.0
openniViewer: /usr/lib/libvtkGraphics.so.5.8.0
openniViewer: /usr/lib/libvtkIO.so.5.8.0
openniViewer: /usr/lib/libvtkFiltering.so.5.8.0
openniViewer: /usr/lib/libvtkCommon.so.5.8.0
openniViewer: /usr/lib/libvtksys.so.5.8.0
openniViewer: CMakeFiles/openniViewer.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable openniViewer"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/openniViewer.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/openniViewer.dir/build: openniViewer
.PHONY : CMakeFiles/openniViewer.dir/build

CMakeFiles/openniViewer.dir/requires: CMakeFiles/openniViewer.dir/src/main.cpp.o.requires
.PHONY : CMakeFiles/openniViewer.dir/requires

CMakeFiles/openniViewer.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/openniViewer.dir/cmake_clean.cmake
.PHONY : CMakeFiles/openniViewer.dir/clean

CMakeFiles/openniViewer.dir/depend:
	cd /home/andrej/studyspace/courseProject/pclTest && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andrej/studyspace/courseProject/pclTest /home/andrej/studyspace/courseProject/pclTest /home/andrej/studyspace/courseProject/pclTest /home/andrej/studyspace/courseProject/pclTest /home/andrej/studyspace/courseProject/pclTest/CMakeFiles/openniViewer.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/openniViewer.dir/depend

