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
  SchemaListBox.h
)

include(FetchContent)

# ===== External Dependency: river
FetchContent_Declare(
  river
  GIT_REPOSITORY    https://github.com/pbotros/river.git
  GIT_TAG           master
  GIT_SHALLOW       1
)
FetchContent_MakeAvailable(river)

# ===== External Dependency: json
FetchContent_Declare(
  json
  GIT_REPOSITORY    https://github.com/nlohmann/json.git
  GIT_TAG           v3.9.1
  GIT_SHALLOW       1
)
FetchContent_MakeAvailable(json)

target_link_libraries(${PLUGIN_NAME} river nlohmann_json::nlohmann_json)

#optional: create IDE groups
plugin_create_filters()
