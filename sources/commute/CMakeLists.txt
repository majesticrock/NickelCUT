add_library(commute_lib STATIC 
    commutator_of_cut.cpp
    restructure_terms.cpp
    extract_flow_coefficients.cpp
)

target_include_directories(commute_lib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(commute_lib 
    PRIVATE 
        mrock::symbolic_operators 
        SymOp_Experimental 
        nickelcut_options
)