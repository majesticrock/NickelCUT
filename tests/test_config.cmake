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

add_executable(momentum_grid_test ${CMAKE_CURRENT_LIST_DIR}/momentum_grid.cpp)
target_link_libraries(momentum_grid_test
    PRIVATE
        shared_lib
        nickelcut_options
)

add_executable(adaptive_flow_test ${CMAKE_CURRENT_LIST_DIR}/adaptive_flow.cpp)
target_compile_definitions(adaptive_flow_test PRIVATE NICKELCUT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(adaptive_flow_test
    PRIVATE
        flow_lib
        nickelcut_options
        mrock::utility
)

add_executable(model_test ${CMAKE_CURRENT_LIST_DIR}/model.cpp)
target_compile_definitions(model_test PRIVATE NICKELCUT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(model_test
    PRIVATE
        shared_lib
        nickelcut_options
        mrock::utility
)


enable_testing()

add_test(NAME momentum_iterator_test COMMAND momentum_iterator_test)
add_test(NAME interaction_data_frame_test COMMAND interaction_data_frame_test)
add_test(NAME momentum_grid_test COMMAND momentum_grid_test)
add_test(NAME adaptive_flow_test COMMAND adaptive_flow_test)
add_test(NAME model_test COMMAND model_test)