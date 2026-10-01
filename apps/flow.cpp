#include "../sources/L.hpp"
#include "../sources/flow/FlowContainer.hpp"
#include "../sources/flow/Model.hpp"
#include "../sources/flow/FlowEquation.hpp"
#include "../sources/helper_functions.hpp"
#include "../sources/flow/numerical_setup.hpp"
#include "../sources/flow/data_file_names.hpp"
#include "../sources/flow/flow_state_serialization.hpp"
#include "../sources/flow/FlowExceptions.hpp"
#include "../sources/flow/ExtractionContainer.hpp"

#define DENSE_BOOK_KEEPING
#ifdef DENSE_BOOK_KEEPING
#include "../sources/flow/DenseBookKeeper.hpp"
typedef NickelCUT::flow::DenseBookKeeper _book_keeper;
#else
#include "../sources/flow/BookKeeper.hpp"
typedef NickelCUT::flow::BookKeeper _book_keeper;
#endif

#include <mrock/utility/InputFileReader.hpp>
#include <mrock/utility/OutputConvenience.hpp>
#include <nlohmann/json.hpp>

#include <cmath>
#include <iostream>
#include <filesystem>

#ifndef OUTPUT_DATA_DIR
#define OUTPUT_DATA_DIR "../../data/nickel_cut/"
#endif

using namespace NickelCUT;
using namespace NickelCUT::flow;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Not enough arguments provided to " << argv[0]
                << "\nUsage: " << argv[0] << "<configfile>" << std::endl;
        return 1;
    }
    mrock::utility::InputFileReader input(argv[1]);

    Model model(input);
    FlowContainer flow_state(model);
    std::cout << "\nConstructed initial states. The filling of the system is " << model.filling << std::endl;

    const std::string output_dir = std::string(OUTPUT_DATA_DIR) 
        + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/") // ensures that OUTPUT_DATA_DIR ends in "/"
        + input.getString("output_dir") + "/"
        + model.data_dir_name();
    const std::string binary_ouput_dir = std::string(OUTPUT_DATA_DIR) 
        + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/") // ensures that OUTPUT_DATA_DIR ends in "/"
        + input.getString("output_dir") + "/"
        + "binaries/"
        + model.data_dir_name();
    
    std::filesystem::create_directories(output_dir);
    std::filesystem::create_directories(binary_ouput_dir);

    FlowEquation flow_equation;
    _book_keeper book_keeper(flow_state, flow_state.band_width(), target_dl(model.U_0), input.getInt("max_runtime"));

    // Default end_reason will probably never be reached.
    std::string end_reason = "Reached l_final";
    double actual_l_final = 0.;
    try {
        boost::numeric::odeint::integrate_adaptive(
                    boost::numeric::odeint::make_controlled<boost_stepper>( abs_error, rel_error ),
                    flow_equation, flow_state, 0.0, l_final, dl(model.U_0), boost::ref(book_keeper));
        actual_l_final = l_final;
    }
    catch (ControlledFlowInterruption& e) {
        actual_l_final = e.get_end_time();
        end_reason = e.what();
        std::cout << end_reason << "\nSaving current state..." << std::endl;
    }

    book_keeper.print_final(flow_state, actual_l_final);

    nlohmann::json j_metadata = model.generate_meta_data_json();
    j_metadata.update(nlohmann::json{{"end_reason", end_reason}});

    nlohmann::json j_flow_data = book_keeper;
    j_flow_data.update(j_metadata);
    nlohmann::json j_final_flow_state = book_keeper.last_good_state;
    j_final_flow_state.update(j_metadata);

#ifdef DENSE_BOOK_KEEPING
    mrock::utility::save_string(j_flow_data.dump(4), output_dir + data_file_names::DENSE_FLOW_STEPS);
#else
    mrock::utility::save_string(j_flow_data.dump(4), output_dir + data_file_names::FLOW_STEPS);
#endif
    
    mrock::utility::save_string(j_final_flow_state.dump(4), output_dir + data_file_names::FINAL_FLOW_STATE_JSON);

    serialize_flow_state(book_keeper.last_good_state, binary_ouput_dir, data_file_names::FINAL_FLOW_STATE_BIN);
    serialize_extracted_channels(ExtractionContainer(book_keeper.last_good_state), binary_ouput_dir, data_file_names::FINAL_EXTRACTED_CHANNELS);

    return 0;
}