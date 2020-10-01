#Open Ephys GUI direcroty-specific file

#add files in this folder
add_sources(open-ephys
	ProcessorGraph.cpp
	ProcessorGraph.h
	ProcessorGraphHttpServer.h
)

#add nested directories
include(FetchContent)

# ===== External Dependency: cpp-httplib
cmake_policy(SET CMP0079 NEW)
set(HTTPLIB_USE_BROTLI_IF_AVAILABLE OFF)
set(HTTPLIB_USE_ZLIB_IF_AVAILABLE OFF)
set(HTTPLIB_USE_OPENSSL_IF_AVAILABLE OFF)
FetchContent_Declare(
  httplib
  GIT_REPOSITORY    https://github.com/yhirose/cpp-httplib.git
  GIT_TAG           v0.7.7
  GIT_SHALLOW       1
)
FetchContent_MakeAvailable(httplib)

# ===== External Dependency: json
FetchContent_Declare(
  json
  GIT_REPOSITORY    https://github.com/nlohmann/json.git
  GIT_TAG           v3.9.1
  GIT_SHALLOW       1
)
FetchContent_MakeAvailable(json)

target_link_libraries(open-ephys httplib::httplib nlohmann_json::nlohmann_json)
