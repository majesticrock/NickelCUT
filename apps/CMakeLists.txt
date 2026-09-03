# NickelCUT/apps/CMakeLists.txt

add_executable(commute commute.cpp)
target_link_libraries(commute 
    PRIVATE 
        nickelcut_options 
        mrock::utility
        mrock::symbolic_operators
        SymOp_Experimental
        commute_lib
)

add_custom_target(run_commute
    COMMAND "${CMAKE_COMMAND}" -E chdir "${CMAKE_SOURCE_DIR}" "$<TARGET_FILE:commute>"
    DEPENDS commute
    USES_TERMINAL
    COMMENT "Running the commute generator"
)

add_executable(flow flow.cpp)
target_link_libraries(flow 
    PRIVATE 
        nickelcut_options
        flow_lib
        mrock::utility
        Boost::serialization
        Boost::iostreams
)

add_executable(resume_flow resume_flow.cpp)
target_link_libraries(resume_flow 
    PRIVATE 
        nickelcut_options
        flow_lib
        mrock::utility
        Boost::serialization
        Boost::iostreams 
)

add_executable(mean_field mean_field.cpp)
target_link_libraries(mean_field
    PRIVATE
        nickelcut_options
        mrock::utility
        mean_field_lib
)