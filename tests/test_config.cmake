include(CTest)

add_executable(wick_test ${CMAKE_CURRENT_LIST_DIR}/wick.cpp)
target_link_libraries(wick_test 
    PRIVATE 
        mrock::symbolic_operators
        nickelcut_options
        SymOp_Experimental
)

add_executable(momentum_iterator_test ${CMAKE_CURRENT_LIST_DIR}/momentum_iterator.cpp)
target_link_libraries(momentum_iterator_test
    PRIVATE
        nickelcut_options
)

add_executable(wick_matrix_test ${CMAKE_CURRENT_LIST_DIR}/wick_matrix.cpp)
target_link_libraries(wick_matrix_test
    PRIVATE
        mrock::symbolic_operators
        nickelcut_options
        SymOp_Experimental
        commute_lib
)

add_executable(wick_second_test ${CMAKE_CURRENT_LIST_DIR}/wick_second.cpp)
target_link_libraries(wick_second_test
    PRIVATE
        mrock::utility
        commute_lib
        nickelcut_options
)


enable_testing()

add_test(NAME wick_test COMMAND wick_test)
add_test(NAME momentum_iterator_test COMMAND momentum_iterator_test)
add_test(NAME wick_matrix_test COMMAND wick_matrix_test)
add_test(NAME wick_second_test COMMAND wick_second_test)