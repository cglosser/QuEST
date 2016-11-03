cmake_minimum_required(VERSION 2.6)
project(quest)
set(CMAKE_CXX_FLAGS "--std=c++11 ${CMAKE_CXX_FLAGS}")

set(quest_VERSION_MAJOR 0)
set(quest_VERSION_MINOR 0)

file(GLOB SOURCES src/*.cpp)

find_package(Eigen3 3.2.2 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Boost 1.55 COMPONENTS program_options REQUIRED)
include_directories(${BOOST_INCLUDE_DIR})

add_executable(quest ${SOURCES})
target_link_libraries(quest ${Boost_LIBRARIES})
