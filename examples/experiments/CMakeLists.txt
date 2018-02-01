add_executable(propagation_error EXCLUDE_FROM_ALL "")
target_sources(propagation_error PUBLIC "${CMAKE_CURRENT_LIST_DIR}/propagation_error.cpp")
target_compile_features(propagation_error PUBLIC ${CPP_FEATURES})
target_link_libraries(propagation_error PUBLIC quest ${Boost_LIBRARIES})

add_executable(comparison EXCLUDE_FROM_ALL "")
target_sources(comparison PUBLIC "${CMAKE_CURRENT_LIST_DIR}/direct_aim_comparison.cpp")
target_compile_features(comparison PUBLIC ${CPP_FEATURES})
target_link_libraries(comparison PUBLIC quest ${Boost_LIBRARIES})

add_executable(aim_direct_SI EXCLUDE_FROM_ALL "")
target_sources(aim_direct_SI PUBLIC "${CMAKE_CURRENT_LIST_DIR}/full_aim_calculation.cpp")
target_compile_features(aim_direct_SI PUBLIC ${CPP_FEATURES})
target_link_libraries(aim_direct_SI PUBLIC quest ${Boost_LIBRARIES})

add_executable(timestepping EXCLUDE_FROM_ALL "")
target_sources(timestepping PUBLIC "${CMAKE_CURRENT_LIST_DIR}/timestepping.cpp")
target_compile_features(timestepping PUBLIC ${CPP_FEATURES})
target_link_libraries(timestepping PUBLIC quest ${Boost_LIBRARIES})
