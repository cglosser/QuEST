set(CMAKE_CXX_FLAGS "--std=c++11 -Wall -fdiagnostics-color=auto ${CMAKE_CXX_FLAGS}")

file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")

find_package(Eigen3 3.2.2 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Boost 1.55 COMPONENTS program_options REQUIRED)
include_directories(${BOOST_INCLUDE_DIR})

add_subdirectory(interactions)

add_executable(quest ${SOURCES})
target_link_libraries(quest ${Boost_LIBRARIES})
