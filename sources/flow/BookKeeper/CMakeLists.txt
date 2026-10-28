add_library(bookkeeper_lib STATIC 
    base.cpp
    Dense.cpp
    Sparse.cpp
    Susceptibility.cpp
)

target_include_directories(bookkeeper_lib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(bookkeeper_lib 
    PRIVATE 
        nickelcut_options
        mrock::utility
    PUBLIC
        shared_lib
        flow_lib
        susceptibility_lib
)