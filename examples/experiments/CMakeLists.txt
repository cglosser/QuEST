add_executable(propagation_error EXCLUDE_FROM_ALL "")
target_sources(propagation_error PUBLIC "${CMAKE_CURRENT_LIST_DIR}/propagation_error.cpp")
target_compile_features(propagation_error PUBLIC ${CPP_FEATURES})
target_link_libraries(propagation_error PUBLIC quest ${Boost_LIBRARIES})

add_executable(aim_check EXCLUDE_FROM_ALL "")
target_sources(aim_check PUBLIC "${CMAKE_CURRENT_LIST_DIR}/aim_check.cpp")
target_compile_features(aim_check PUBLIC ${CPP_FEATURES})
target_link_libraries(aim_check PUBLIC quest ${Boost_LIBRARIES})

add_executable(multiphysics EXCLUDE_FROM_ALL "")
target_sources(multiphysics PUBLIC "${CMAKE_CURRENT_LIST_DIR}/multiphysics.cpp")
target_compile_features(multiphysics PUBLIC ${CPP_FEATURES})
target_link_libraries(multiphysics PUBLIC quest ${Boost_LIBRARIES})

add_executable(timestepping EXCLUDE_FROM_ALL "")
target_sources(timestepping PUBLIC "${CMAKE_CURRENT_LIST_DIR}/timestepping.cpp")
target_compile_features(timestepping PUBLIC ${CPP_FEATURES})
target_link_libraries(timestepping PUBLIC quest ${Boost_LIBRARIES})

add_executable(hotfix EXCLUDE_FROM_ALL "")
target_sources(hotfix PUBLIC "${CMAKE_CURRENT_LIST_DIR}/hotfix.cpp")
target_compile_features(hotfix PUBLIC ${CPP_FEATURES})
target_link_libraries(hotfix PUBLIC quest ${Boost_LIBRARIES})

add_executable(hotfix_SI EXCLUDE_FROM_ALL "")
target_sources(hotfix_SI PUBLIC "${CMAKE_CURRENT_LIST_DIR}/hotfix_SI.cpp")
target_compile_features(hotfix_SI PUBLIC ${CPP_FEATURES})
target_link_libraries(hotfix_SI PUBLIC quest ${Boost_LIBRARIES})

add_executable(crossover EXCLUDE_FROM_ALL "")
target_sources(crossover PUBLIC "${CMAKE_CURRENT_LIST_DIR}/crossover.cpp")
target_compile_features(crossover PUBLIC ${CPP_FEATURES})
target_link_libraries(crossover PUBLIC quest ${Boost_LIBRARIES})

add_executable(potential_and_field_error EXCLUDE_FROM_ALL "")
target_sources(potential_and_field_error PUBLIC "${CMAKE_CURRENT_LIST_DIR}/potential_and_field_error.cpp")
target_compile_features(potential_and_field_error PUBLIC ${CPP_FEATURES})
target_link_libraries(potential_and_field_error PUBLIC quest ${Boost_LIBRARIES})
