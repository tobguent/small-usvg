set(EXECUTABLE_NAME demo)

add_definitions(-DNOMINMAX)
add_executable(${EXECUTABLE_NAME} main.cpp)
target_link_libraries(${EXECUTABLE_NAME} PRIVATE usvg)

target_compile_definitions(demo PRIVATE "USVG_SCENES_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/scenes/\"")

# IDE Support
# --------------------------------------------------
if (USVG_IDE_SUPPORT)
   enable_ide_support_vislab(${EXECUTABLE_NAME} "" ${CMAKE_CURRENT_SOURCE_DIR})
endif()