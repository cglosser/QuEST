include("${CMAKE_CURRENT_LIST_DIR}/integrator/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/interactions/CMakeLists.txt")

target_sources(quest
  PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/common.h"
    "${CMAKE_CURRENT_LIST_DIR}/configuration.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/configuration.h"
    "${CMAKE_CURRENT_LIST_DIR}/fourier.h"
    "${CMAKE_CURRENT_LIST_DIR}/lagrange_set.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/lagrange_set.h"
    "${CMAKE_CURRENT_LIST_DIR}/math_utils.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/math_utils.h"
    "${CMAKE_CURRENT_LIST_DIR}/pulse.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/pulse.h"
    "${CMAKE_CURRENT_LIST_DIR}/quantum_dot.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/quantum_dot.h"
)

add_executable(quest-bin "")
set_target_properties(quest-bin PROPERTIES OUTPUT_NAME quest)

target_sources(quest-bin
  PUBLIC
    "${CMAKE_CURRENT_LIST_DIR}/main.cpp"
)

target_compile_features(quest-bin PRIVATE ${CPP_FEATURES})
target_link_libraries(quest-bin PUBLIC quest ${Boost_LIBRARIES})
