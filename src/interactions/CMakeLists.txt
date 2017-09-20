include("${CMAKE_CURRENT_LIST_DIR}/AIM/CMakeLists.txt")

target_sources(quest
  PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/direct_interaction.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/direct_interaction.h"
    "${CMAKE_CURRENT_LIST_DIR}/green_function.h"
    "${CMAKE_CURRENT_LIST_DIR}/history_interaction.h"
    "${CMAKE_CURRENT_LIST_DIR}/interaction.h"
    "${CMAKE_CURRENT_LIST_DIR}/pulse_interaction.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/pulse_interaction.h"
)
