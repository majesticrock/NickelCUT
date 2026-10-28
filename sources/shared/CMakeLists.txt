add_library(shared_lib STATIC 
    Model.cpp
    occupation_numbers.cpp
)

target_include_directories(shared_lib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(shared_lib 
    PRIVATE 
        nickelcut_options
        mrock::utility
)