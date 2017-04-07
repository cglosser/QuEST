find_package(Eigen3 3.2.2 REQUIRED)
target_include_directories(quest PUBLIC ${EIGEN3_INCLUDE_DIR})

find_package(Boost 1.55 COMPONENTS program_options REQUIRED)
target_include_directories(quest PUBLIC ${Boost_INCLUDE_DIR})
target_link_libraries(quest PUBLIC ${Boost_LIBRARIES})

include(src/interactions/CMakeLists.txt)

target_sources(quest
  PUBLIC
    "${CMAKE_CURRENT_LIST_DIR}/common.h"
    "${CMAKE_CURRENT_LIST_DIR}/configuration.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/configuration.h"
    "${CMAKE_CURRENT_LIST_DIR}/history.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/history.h"
    "${CMAKE_CURRENT_LIST_DIR}/integrator.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/integrator.h"
    "${CMAKE_CURRENT_LIST_DIR}/lagrange_set.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/lagrange_set.h"
    "${CMAKE_CURRENT_LIST_DIR}/main.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/math_utils.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/math_utils.h"
    "${CMAKE_CURRENT_LIST_DIR}/pulse.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/pulse.h"
    "${CMAKE_CURRENT_LIST_DIR}/quantum_dot.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/quantum_dot.h"
)
