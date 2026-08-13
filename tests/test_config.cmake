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

enable_testing()

add_test(NAME wick_test COMMAND wick_test)
add_test(NAME momentum_iterator_test COMMAND momentum_iterator_test)