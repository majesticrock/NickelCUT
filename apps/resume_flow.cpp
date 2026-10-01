#include "../sources/L.hpp"
#include "../sources/flow/FlowContainer.hpp"
#include "../sources/flow/Model.hpp"
#include "../sources/flow/FlowEquation.hpp"
#include "../sources/helper_functions.hpp"
#include "../sources/flow/BookKeeper.hpp"
#include "../sources/flow/numerical_setup.hpp"
#include "../sources/flow/data_file_names.hpp"
#include "../sources/flow/flow_state_serialization.hpp"
#include "../sources/flow/FlowExceptions.hpp"
#include "../sources/flow/ExtractionContainer.hpp"

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
    if (argc < 3) {
        std::cerr << "Not enough arguments provided to " << argv[0]
                << "\nUsage: " << argv[0] << "<configfile> <int: resume_step>" << std::endl;
        return 1;
    }

    const bool resume_run = argc >= 3;
    const int resume_step = resume_run ? std::stoi(argv[2]) : 0;
    (void)resume_step;

    mrock::utility::InputFileReader input(argv[1]);
    Model model(input);

    const std::string output_dir = std::string(OUTPUT_DATA_DIR) 
        + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/")
        + input.getString("output_dir") + "/"
        + model.data_dir_name();
    const std::string binary_ouput_dir = std::string(OUTPUT_DATA_DIR) 
        + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/")
        + input.getString("output_dir") + "/"
        + "binaries/"
        + model.data_dir_name();

    std::filesystem::create_directories(output_dir);
    std::filesystem::create_directories(binary_ouput_dir);

    FlowContainer flow_state;
    if (resume_run) {
        const std::string final_state_file = binary_ouput_dir + data_file_names::FINAL_FLOW_STATE;
        if (!std::filesystem::exists(final_state_file)) {
            std::cerr << "No serialized flow state found at " << final_state_file << ". Cannot resume flow." << std::endl;
            return 1;
        }
        flow_state = deserialize_flow_state(binary_ouput_dir, data_file_names::FINAL_FLOW_STATE);
        std::cout << "Loaded serialized flow state from " << final_state_file << std::endl;
    }
    else {
        flow_state = FlowContainer(model);
        std::cout << "\nConstructed initial states. The filling of the system is " << model.filling << std::endl;
    }

    FlowEquation flow_equation;
    BookKeeper book_keeper(flow_state, flow_state.band_width(), target_dl(model.U_0), input.getInt("max_runtime"));

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

    if (!(flow_state.is_inversion_symmetric() && flow_state.is_hermitian())) {
        std::cerr << "State is no longer reliable!" << std::endl;
    }

    book_keeper.print_final(flow_state, actual_l_final);

    nlohmann::json j_metadata = model.generate_meta_data_json();
    j_metadata.update(nlohmann::json{{"end_reason", end_reason}});

    nlohmann::json j_flow_data = book_keeper;
    j_flow_data.update(j_metadata);
    j_flow_data.update({"extracted_channels", ExtractionContainer(flow_state)});

    nlohmann::json j_full_flow_state = flow_state;
    j_full_flow_state.update(j_metadata);

    mrock::utility::save_string(j_flow_data.dump(4), output_dir + data_file_names::FLOW_STEPS);
    mrock::utility::save_string(j_full_flow_state.dump(4), output_dir + data_file_names::FULL_FLOW_STATE);
    serialize_flow_state(flow_state, binary_ouput_dir, data_file_names::FINAL_FLOW_STATE);
    serialize_extracted_channels(ExtractionContainer(flow_state), binary_ouput_dir, data_file_names::EXTRACTED_CHANNELS);

    return 0;
}