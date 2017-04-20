find_package (Boost COMPONENTS unit_test_framework)
if(${Boost_FOUND})
  add_executable(qtest "")
  target_compile_features(qtest PUBLIC ${CPP_FEATURES})
  target_sources(qtest PRIVATE "${CMAKE_CURRENT_LIST_DIR}/main.cpp")
  target_sources(qtest
    PUBLIC
      "${CMAKE_CURRENT_LIST_DIR}/history_test.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/lagrange_set_test.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/quantum_dot_test.cpp"
  )
  target_link_libraries(qtest PUBLIC libquest ${Boost_LIBRARIES})
endif(${Boost_FOUND})
