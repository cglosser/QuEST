cmake_minimum_required(VERSION 3.1)
project(QuEST VERSION 0.1 LANGUAGES CXX)

# Location of custom detection scripts
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

# Set RELEASE as default build type
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Enable warning flags if available
include("${PROJECT_SOURCE_DIR}/cmake/flag_check.cmake")
enable_cxx_compiler_flag_if_supported("-Wall")
enable_cxx_compiler_flag_if_supported("-Wextra")
enable_cxx_compiler_flag_if_supported("-pedantic")

add_executable(quest "")
target_compile_features(quest PRIVATE
  cxx_auto_type
  cxx_constexpr
  cxx_delegating_constructors
  cxx_lambdas
  cxx_range_for
  cxx_strong_enums
)

include(src/CMakeLists.txt)
