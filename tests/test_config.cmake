include(CTest)

add_executable(wick_test ${CMAKE_CURRENT_LIST_DIR}/wick.cpp)
target_link_libraries(wick_test 
    PRIVATE 
        mrock::symbolic_operators
        nickelcut_options
        SymOp_Experimental
)

enable_testing()

add_test(NAME wick_test COMMAND wick_test)