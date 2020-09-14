#plugin build file
cmake_minimum_required(VERSION 3.5.0)

#include common rules
include(../PluginRules.cmake)

#add sources, not including OpenEphysLib.cpp
add_sources(${PLUGIN_NAME}
	SpikeSorter.cpp
	SpikeSorter.h
	SpikeSortBoxes.cpp
	SpikeSortBoxes.h
	SpikeSorterEditor.cpp
	SpikeSorterEditor.h
	SpikeSorterCanvas.cpp
	SpikeSorterCanvas.h
	)


# ===== External Dependency: TCrosser
configure_file(CMakeLists.txt.in.tcrosser tcrosser-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
		RESULT_VARIABLE result
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tcrosser-download )
if(result)
	message(FATAL_ERROR "CMake step for tcrosser failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
		RESULT_VARIABLE result
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tcrosser-download )
if(result)
	message(FATAL_ERROR "Build step for tcrosser failed: ${result}")
endif()
set(BUILD_TCROSSER ON)
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/tcrosser-src
		${CMAKE_CURRENT_BINARY_DIR}/tcrosser-build
		EXCLUDE_FROM_ALL)

#optional: create IDE groups
#plugin_create_filters()

target_link_libraries(${PLUGIN_NAME} tcrosser)
