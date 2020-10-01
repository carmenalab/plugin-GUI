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


include(FetchContent)

# ===== External Dependency: TCrosser
FetchContent_Declare(
  tcrosser
  GIT_REPOSITORY    git@github.com:carmenalab/c.git
  GIT_TAG           master
  GIT_SHALLOW       1
)
FetchContent_MakeAvailable(tcrosser)

# ===== External Dependency: json
FetchContent_Declare(
  json
  GIT_REPOSITORY    https://github.com/nlohmann/json.git
  GIT_TAG           v3.9.1
  GIT_SHALLOW       1
)
FetchContent_MakeAvailable(json)

#optional: create IDE groups
#plugin_create_filters()

target_link_libraries(${PLUGIN_NAME} tcrosser nlohmann_json::nlohmann_json)
