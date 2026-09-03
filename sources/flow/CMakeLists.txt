option(NICKELCUT_RUN_COMMUTE "Build and run the commute target to regenerate FlowEquation.cpp" OFF)
option(NICKELCUT_FORCE_FLOW_EQUATION_REGEN "Always regenerate FlowEquation.cpp before building flow_lib" OFF)

set(FLOW_EQUATION_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/FlowEquation.cpp")

if(NICKELCUT_RUN_COMMUTE OR NICKELCUT_FORCE_FLOW_EQUATION_REGEN OR NOT EXISTS "${FLOW_EQUATION_SOURCE}")
    add_custom_command(
        OUTPUT "${FLOW_EQUATION_SOURCE}"
        COMMAND "${CMAKE_COMMAND}" -E chdir "${CMAKE_SOURCE_DIR}" "$<TARGET_FILE:commute>"
        DEPENDS commute
        VERBATIM
    )
    add_custom_target(generate_flow_equation ALL DEPENDS "${FLOW_EQUATION_SOURCE}")
endif()

add_library(flow_lib STATIC
    BookKeeper.cpp
    DecouplingChannel.cpp
    Model.cpp
    FlowContainer.cpp
    occupation_numbers.cpp
    flow_state_serialization.cpp
    "${FLOW_EQUATION_SOURCE}"
)

if(TARGET generate_flow_equation)
    add_dependencies(flow_lib generate_flow_equation)
endif()

target_include_directories(flow_lib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(flow_lib 
    PRIVATE 
        nickelcut_options
        mrock::utility
    PUBLIC
        OpenMP::OpenMP_CXX 
)