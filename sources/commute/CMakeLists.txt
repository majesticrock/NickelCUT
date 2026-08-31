add_library(commute_lib STATIC 
    commutator_of_cut.cpp
    
    Hamiltonian/restructure_terms.cpp
    Hamiltonian/extract_flow_coefficients.cpp
    Hamiltonian/export_as_flow_equation.cpp
    Hamiltonian/commute_and_extract.cpp

    PairAnnihilation/restructure_terms.cpp
    PairAnnihilation/commute_and_extract.cpp

    verify/Verifier.cpp
)

target_include_directories(commute_lib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(commute_lib 
    PRIVATE  
        nickelcut_options
    PUBLIC
        SymOp_Experimental
        mrock::symbolic_operators 
)