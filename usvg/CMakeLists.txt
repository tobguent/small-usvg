# CMake Library Setup
# --------------------------------------------------
set(LIBRARY_NAME usvg)

# Find source files
file(GLOB SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

# Create library
add_library(${LIBRARY_NAME} ${SOURCES})

# link
target_link_libraries(${LIBRARY_NAME} PUBLIC eigen3::eigen tinyxml2)

# include directories
# public:
target_include_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# private:
target_include_directories(${LIBRARY_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# External Dependency Setup
# --------------------------------------------------
include(cmake/eigen3.cmake)
include(cmake/openmp.cmake)
include(cmake/tinyxml2.cmake)

# IDE Support
# --------------------------------------------------
if (USVG_IDE_SUPPORT)
   enable_ide_support_vislab(${LIBRARY_NAME} "" ${CMAKE_CURRENT_SOURCE_DIR})
endif()
