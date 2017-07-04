include("${CMAKE_CURRENT_LIST_DIR}/RHS/CMakeLists.txt")

target_sources(quest
  PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/history.h"
    "${CMAKE_CURRENT_LIST_DIR}/integrator.h"
    "${CMAKE_CURRENT_LIST_DIR}/weights.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/weights.h"
)
