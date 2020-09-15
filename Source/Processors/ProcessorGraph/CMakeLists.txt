#Open Ephys GUI direcroty-specific file

#add files in this folder
add_sources(open-ephys
	ProcessorGraph.cpp
	ProcessorGraph.h
	ProcessorGraphHttpServer.h
)

#add nested directories

# ===== External Dependency: cpp-httplib
cmake_policy(SET CMP0079 NEW)
set(HTTPLIB_USE_BROTLI_IF_AVAILABLE OFF)
set(HTTPLIB_USE_ZLIB_IF_AVAILABLE OFF)

configure_file(CMakeLists.txt.in.httplib httplib-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
		RESULT_VARIABLE result
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/httplib-download )
if(result)
	message(FATAL_ERROR "CMake step for httplib failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
		RESULT_VARIABLE result
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/httplib-download )
if(result)
	message(FATAL_ERROR "Build step for httplib failed: ${result}")
endif()
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/httplib-src
		${CMAKE_CURRENT_BINARY_DIR}/httplib-build
		EXCLUDE_FROM_ALL)

# ===== External Dependency: json
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


target_link_libraries(open-ephys httplib::httplib nlohmann_json::nlohmann_json)
