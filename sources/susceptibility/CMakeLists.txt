add_library(susceptibility_lib STATIC 
    ieom_base.cpp
    SC.cpp
    Spin.cpp
    Charge.cpp
)

target_include_directories(susceptibility_lib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(susceptibility_lib 
    PRIVATE 
        nickelcut_options
        mrock::utility
    PUBLIC
        shared_lib
)