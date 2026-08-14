add_library(mean_field_lib STATIC 
    Model.cpp
)

target_include_directories(mean_field_lib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(mean_field_lib 
    PRIVATE 
        nickelcut_options
        Boost::serialization
        Boost::iostreams
)