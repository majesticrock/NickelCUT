#include "../sources/flow/FlowContainer.hpp"
#include "../sources/shared/Model.hpp"
#include "../sources/flow/data_file_names.hpp"
#include "../sources/flow/flow_state_serialization.hpp"
#include "../sources/susceptibility/Charge.hpp"
#include "../sources/susceptibility/SC.hpp"
#include "../sources/susceptibility/Spin.hpp"

#include <mrock/utility/InputFileReader.hpp>
#include <mrock/utility/OutputConvenience.hpp>
#include <nlohmann/json.hpp>

#include <iostream>

#ifndef OUTPUT_DATA_DIR
#define OUTPUT_DATA_DIR "build/"  //"../../data/nickel_cut/"
#endif

using namespace NickelCUT;
using namespace NickelCUT::susceptibility;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Invalid number of arguments: Use <path_to_executable> <configfile>" << std::endl;
        return -1;
    }
    mrock::utility::InputFileReader input(argv[1]);
    const Model model(input);

    const std::string load_binary_dir =
        std::string(OUTPUT_DATA_DIR) +
        (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/")  // ensures that OUTPUT_DATA_DIR ends in "/"
        + input.getString("output_dir") + "/" + "binaries/";
    const flow::FlowContainer final_flow_state = flow::deserialize_flow_state(
        load_binary_dir + model.data_dir_name(), flow::data_file_names::FINAL_FLOW_STATE_BIN);

    const std::string output_dir =
        std::string(OUTPUT_DATA_DIR) +
        (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/")  // ensures that OUTPUT_DATA_DIR ends in "/"
        + input.getString("output_dir") + "/" + model.data_dir_name();

    const nlohmann::json j_metadata = model.generate_meta_data_json();

    SC ieom_sc(final_flow_state, model);
    Spin ieom_spin(&ieom_sc);
    Charge ieom_charge(&ieom_sc);

    const std::vector<ieom_base::Vector> states = {ieom_sc.swave_state(), ieom_sc.dwave_state()};

    const auto sc_result = ieom_sc.static_susceptibility_at_high_symmetry(states);
    const auto spin_result = ieom_spin.static_susceptibility_at_high_symmetry(states);
    const auto charge_result = ieom_charge.static_susceptibility_at_high_symmetry(states);

    nlohmann::json jData = {
        {"SC", sc_result},
        {"spin", spin_result},
        {"charge", charge_result},
        {"dispersion", ieom_sc.get_dispersion()},
        {"negative_excitation_sc", ieom_sc.found_negative_excitation()},
        {"negative_excitation_spin", ieom_spin.found_negative_excitation()},
        {"negative_excitation_charge", ieom_charge.found_negative_excitation()},
        {"pomeranchuk_instability", ieom_sc.found_pomeranchuk_instability()},
        {"hartree_fock_filling", ieom_sc.hartree_fock_filling()}
    };
    jData.update(j_metadata);

    mrock::utility::save_string(jData.dump(4), output_dir + "ieom_susceptibility.json.gz");
    
    return 0;
}