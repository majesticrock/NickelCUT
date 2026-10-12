include(CTest)

add_executable(momentum_iterator_test ${CMAKE_CURRENT_LIST_DIR}/momentum_iterator.cpp)
target_link_libraries(momentum_iterator_test
    PRIVATE
        flow_lib
        nickelcut_options
)

add_executable(interaction_data_frame_test ${CMAKE_CURRENT_LIST_DIR}/interaction_data_frame.cpp)
target_link_libraries(interaction_data_frame_test
    PRIVATE
        flow_lib
        nickelcut_options
)


enable_testing()

add_test(NAME momentum_iterator_test COMMAND momentum_iterator_test)
add_test(NAME interaction_data_frame_test COMMAND interaction_data_frame_test)