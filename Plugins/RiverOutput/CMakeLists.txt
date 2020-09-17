#plugin build file
cmake_minimum_required(VERSION 3.5.0)
set(CMAKE_CXX_STANDARD 14)

#include common rules
include(../PluginRules.cmake)

#add sources, not including OpenEphysLib.cpp
add_sources(${PLUGIN_NAME}
  RiverOutput.cpp
  RiverOutput.h
  RiverOutputEditor.cpp
  RiverOutputEditor.h
  )

# ===== External Dependency: river
configure_file(CMakeLists.txt.in.river river-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/river-download )
if(result)
  message(FATAL_ERROR "CMake step for river failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/river-download )
if(result)
  message(FATAL_ERROR "Build step for river failed: ${result}")
endif()
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/river-src
  ${CMAKE_CURRENT_BINARY_DIR}/river-build
  EXCLUDE_FROM_ALL)

# ===== External Dependency: json
if (NOT TARGET nlohmann_json::nlohmann_json)
  configure_file(CMakeLists.txt.in.json json-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/json-download )
  if(result)
    message(FATAL_ERROR "CMake step for json failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/json-download )
  if(result)
    message(FATAL_ERROR "Build step for json failed: ${result}")
  endif()
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/json-src
          ${CMAKE_CURRENT_BINARY_DIR}/json-build
          EXCLUDE_FROM_ALL)
endif()

target_link_libraries(${PLUGIN_NAME} river nlohmann_json::nlohmann_json)

#optional: create IDE groups
plugin_create_filters()
