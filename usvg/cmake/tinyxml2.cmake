FetchContent_Declare(
	tinyxml2
	GIT_REPOSITORY https://github.com/leethomason/tinyxml2.git
)
set(tinyxml2_BUILD_TESTING OFF CACHE BOOL "Build tests for tinyxml2")
FetchContent_MakeAvailable(tinyxml2)
set_target_properties(
	tinyxml2
	PROPERTIES
	ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
	LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_target_properties(tinyxml2 PROPERTIES FOLDER "extern")
