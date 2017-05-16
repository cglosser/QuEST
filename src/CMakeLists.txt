include("${CMAKE_CURRENT_LIST_DIR}/interactions/CMakeLists.txt")

target_sources(libquest
  PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/common.h"
    "${CMAKE_CURRENT_LIST_DIR}/configuration.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/configuration.h"
    "${CMAKE_CURRENT_LIST_DIR}/history.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/history.h"
    "${CMAKE_CURRENT_LIST_DIR}/integrator.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/integrator.h"
    "${CMAKE_CURRENT_LIST_DIR}/lagrange_set.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/lagrange_set.h"
    "${CMAKE_CURRENT_LIST_DIR}/math_utils.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/math_utils.h"
    "${CMAKE_CURRENT_LIST_DIR}/pulse.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/pulse.h"
    "${CMAKE_CURRENT_LIST_DIR}/quantum_dot.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/quantum_dot.h"
)

add_executable(quest "")

target_sources(quest
  PUBLIC
    "${CMAKE_CURRENT_LIST_DIR}/main.cpp"
)

target_compile_features(quest PRIVATE ${CPP_FEATURES})
target_link_libraries(quest PUBLIC libquest ${Boost_LIBRARIES})
