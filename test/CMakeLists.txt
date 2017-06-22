find_package (Boost COMPONENTS unit_test_framework)
if(${Boost_FOUND})
  add_executable(qtest-bin "")
  set_target_properties(qtest-bin PROPERTIES OUTPUT_NAME qtest)

  target_compile_features(qtest-bin PUBLIC ${CPP_FEATURES})
  target_sources(qtest-bin PRIVATE "${CMAKE_CURRENT_LIST_DIR}/main.cpp")
  target_sources(qtest-bin
    PUBLIC
      "${CMAKE_CURRENT_LIST_DIR}/integrator_test.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/lagrange_set_test.cpp"
<<<<<<< HEAD
      "${CMAKE_CURRENT_LIST_DIR}/quantum_dot_test.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/history_interaction_test.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/pulse_interaction_test.cpp"
      )
  target_link_libraries(qtest PUBLIC libquest ${Boost_LIBRARIES})
=======
      "${CMAKE_CURRENT_LIST_DIR}/pulse_interaction_test.cpp"
      "${CMAKE_CURRENT_LIST_DIR}/quantum_dot_test.cpp"
  )
  target_link_libraries(qtest-bin PUBLIC quest ${Boost_LIBRARIES})
>>>>>>> 9ed9f70a7fc64e18c1d413870ab551425712af26
endif(${Boost_FOUND})
