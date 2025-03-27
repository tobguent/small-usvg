cmake_minimum_required(VERSION 3.18)

project(small-usvg
        LANGUAGES CXX
        DESCRIPTION "Example implementation of the paper Unified Smooth Vector Graphics: Modeling Gradient Meshes and Curve-based Approaches Jointly as Poisson Problem."
        HOMEPAGE_URL "https://vc.tf.fau.de/")

# ------------------------------------
# ------- Helper Functions -----------
# ------------------------------------
include(FetchContent)
include(cmake/enableIdeSupport.cmake)

# ------------------------------------
# ------------ Options ---------------
# ------------------------------------

# The following options' default value depends on whether we are using a multi config generator such as Visual Studio or
# or XCode. We first have to load that global property into a local variable to be able to read it.
get_property(IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

if (IS_MULTI_CONFIG)
    option(USVG_IDE_SUPPORT
            "Enable support for IDEs such as Visual Studio or XCode (targets grouped into folders and source groups specified)." ON)
else()
    option(USVG_IDE_SUPPORT
            "Enable support for IDEs such as Visual Studio or XCode (targets grouped into folders and source groups specified)." OFF)
endif()

if (IS_MULTI_CONFIG OR CMAKE_BUILD_TYPE MATCHES Release)
    option(USVG_USE_OPENMP "Enable OpenMP support?" ON)
else()
    option(USVG_USE_OPENMP "Enable OpenMP support?" OFF)
endif()

# ------------------------------------
# ----- Compiler settings ------------
# ------------------------------------
# C++17
set(CMAKE_CXX_STANDARD 17)

if(WIN32)
  add_definitions(-DNOMINMAX)
  add_compile_options(-bigobj)
endif()

# ------------------------------------
# ------- CMake Settings -------------
# ------------------------------------

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (USVG_IDE_SUPPORT)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

# ------------------------------------
# ----- Compile usvg library ---------
# ------------------------------------
add_subdirectory(usvg)

# ------------------------------------
# ----- Compile demo app -------------
# ------------------------------------
add_subdirectory(demo)
