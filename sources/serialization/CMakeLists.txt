add_library(serialization_lib STATIC 
    ExtractionContainer.cpp
    FlowContainer.cpp
)

target_include_directories(serialization_lib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(serialization_lib 
    PRIVATE 
        nickelcut_options
        flow_lib
        shared_lib
        Boost::serialization
        Boost::iostreams
)