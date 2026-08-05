add_library(SymOp_Experimental STATIC 
    WickOrderedCollector.cpp
    WickOrderedExpression.cpp
    WickOrderedTerm.cpp
)

target_include_directories(SymOp_Experimental PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(SymOp_Experimental PRIVATE mrock::symbolic_operators nickelcut_options)