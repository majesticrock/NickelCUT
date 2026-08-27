cmake_minimum_required(VERSION 3.30)

project(NickelCUT VERSION 2.1.0 LANGUAGES CXX)

# Keep Release optimizations while leaving assertions enabled.
set(CMAKE_CXX_FLAGS_BETA_TEST "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING
    "Flags for an optimized build with assertions enabled." FORCE)
string(REGEX REPLACE "(^| )[/-]DNDEBUG( |$)" "\\1\\2"
    CMAKE_CXX_FLAGS_BETA_TEST "${CMAKE_CXX_FLAGS_BETA_TEST}")

# Set local user-specific settings
# By default, this adds ~/usr/local/ and ../../.mrock/ to the cmake and include paths
# If additional configurations/includes are needed, provide them in cmake/user_settings.cmake
# Extra include dirs are stored in the cmake variable EXTRA_INCLUDE_DIRS
include(cmake/include_settings.cmake)

# Build type settings; defaults to Release
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo" "BetaTest" "MinSizeRel")
endif()
message(STATUS "Build type is '${CMAKE_BUILD_TYPE}'.")



#
# Library specifics
#

find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(
    Boost REQUIRED
    COMPONENTS
        serialization
        iostreams
)
find_package(mrock 
    REQUIRED 
    COMPONENTS 
        symbolic_operators 
        utility
)

# provided by the mrock library
mrock_generate_information_header(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/nickelcut/info.h"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUT_HEADER MROCK_GENERATED_INFO_HEADER
    OUT_INCLUDE_DIR MROCK_GENERATED_INCLUDE_DIR
)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/nickelcut_options.cmake)
# Sets include dir for the generated metadata header file
target_include_directories(nickelcut_options INTERFACE ${MROCK_GENERATED_INCLUDE_DIR})

add_subdirectory(sources/experimental)
add_subdirectory(sources/commute)
add_subdirectory(sources/flow)
add_subdirectory(sources/mean_field)

add_subdirectory(apps)

option(BUILD_TESTING 
    "Build and configure the testing suit" 
    OFF
)

if(BUILD_TESTING)
    include(tests/test_config.cmake)
endif()