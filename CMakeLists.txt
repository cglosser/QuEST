cmake_minimum_required(VERSION 2.6)
project(quest)

# Location of custom detection scripts
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

set(CMAKE_CXX_FLAGS "--std=c++11 -Wall -fdiagnostics-color=auto ${CMAKE_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")

set(quest_VERSION_MAJOR 0)
set(quest_VERSION_MINOR 0)

find_package(Eigen3 3.2.2 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Boost 1.55 COMPONENTS program_options REQUIRED)
include_directories(${BOOST_INCLUDE_DIR})

add_subdirectory(src)
add_subdirectory(src/interactions)

add_executable(quest ${SOURCES})
target_link_libraries(quest ${Boost_LIBRARIES})
