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
    mrock::utility::InputFileReader input(argv[1]);
    const int resume_step = std::stoi(argv[2]);
    Model model(input);

    const std::string output_dir = std::string(OUTPUT_DATA_DIR) 
        + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/") // ensures that OUTPUT_DATA_DIR ends in "/"
        + input.getString("output_dir") + "/"
        + model.data_dir_name();
    const std::string binary_ouput_dir = std::string(OUTPUT_DATA_DIR) 
        + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/") // ensures that OUTPUT_DATA_DIR ends in "/"
        + input.getString("output_dir") + "/"
        + "binaries/"
        + model.data_dir_name();

    if (!(std::filesystem::exists(output_dir) && std::filesystem::exists(binary_ouput_dir))) {
        std::cerr << "Output dir does not exist. Cannot resume flow." << std::endl;
    }

    FlowContainer flow_state = deserialize_flow_state(binary_ouput_dir, data_file_names::FINAL_FLOW_STATE + (resume_step > 0 ? argv[1] : ""));
    FlowEquation flow_equation;
    BookKeeper book_keeper(flow_state, target_dl(model.U_0), input.getInt("max_runtime"));

    double actual_l_final = -1.;
    try {
        // Since the flow equation does not explicitly depend on l, we can just tell it to start again at l=0 and go to l_final
        // while keeping in mind that l=0 now corresponds to the l at which the last computation ended.
        boost::numeric::odeint::integrate_adaptive(
                    boost::numeric::odeint::make_controlled<boost_stepper>( abs_error, rel_error ),
                    flow_equation, flow_state, 0.0, l_final(model.U_0), dl(model.U_0), boost::ref(book_keeper));
        actual_l_final = l_final(model.U_0);
    }
    catch (ControlledFlowInterruption& e) {
        actual_l_final = e.get_end_time();
        std::cout << e.what() << "\nSaving current state..." << std::endl;
    }

    // Checks whether symmetries are preserved
    if (!(book_keeper.lowest_ROD_state.is_inversion_symmetric() && book_keeper.lowest_ROD_state.is_hermitian())) {
        std::cerr << "State is no longer reliable!" << std::endl;
    }

    book_keeper.print_final(flow_state, actual_l_final);

    const nlohmann::json j_metadata = model.generate_meta_data_json();
    nlohmann::json j_flow_data = book_keeper;
    j_flow_data.merge_patch(j_metadata);
    nlohmann::json j_full_flow_state = book_keeper.lowest_ROD_state;
    j_full_flow_state.merge_patch(j_metadata);

    const std::string name_append = std::to_string(resume_step + 1);
    mrock::utility::save_string(j_flow_data.dump(4), output_dir + data_file_names::FLOW_STEPS + name_append);
    mrock::utility::save_string(j_full_flow_state.dump(4), output_dir + data_file_names::FULL_FLOW_STATE + name_append);
    serialize_flow_state(book_keeper.lowest_ROD_state, binary_ouput_dir, data_file_names::LOWEST_ROD_STATE + name_append);
    serialize_flow_state(flow_state, binary_ouput_dir, data_file_names::FINAL_FLOW_STATE + name_append);
    serialize_extracted_channels(book_keeper.extracted_channels[book_keeper.index_of_lowest_ROD], binary_ouput_dir, data_file_names::LOWEST_ROD_EXTRACTED_CHANNELS + name_append);

    return 0;
}