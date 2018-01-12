add_executable(propagation_error EXCLUDE_FROM_ALL "")
target_sources(propagation_error PUBLIC "${CMAKE_CURRENT_LIST_DIR}/propagation_error.cpp")
target_compile_features(propagation_error PUBLIC ${CPP_FEATURES})
target_link_libraries(propagation_error PUBLIC quest ${Boost_LIBRARIES})
