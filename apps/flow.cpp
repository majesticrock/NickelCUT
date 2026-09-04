#include "../sources/L.hpp"
#include "../sources/flow/momentum_iterator.hpp"
#include "../sources/flow/FlowContainer.hpp"
#include "../sources/flow/Model.hpp"
#include "../sources/flow/FlowEquation.hpp"
#include "../sources/helper_functions.hpp"
#include "../sources/flow/BookKeeper.hpp"
#include "../sources/flow/numerical_setup.hpp"
#include "../sources/flow/data_file_names.hpp"
#include "../sources/flow/flow_state_serialization.hpp"

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

int main(int /*argc*/, char** /*argv*/) {
    Model model(U_0, tprime, mu_0, T);
    FlowContainer flow_state(model);
    std::cout << "\nConstructed initial states. The filling of the system is " << model.filling << std::endl;

    const std::string output_folder = std::string(OUTPUT_DATA_DIR) 
        + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/") // ensures that OUTPUT_DATA_DIR ends in "/"
        + model.data_dir_name();
    std::filesystem::create_directories(output_folder);

    FlowEquation flow_equation;
    BookKeeper book_keeper(flow_state, target_dl);

    try {
        boost::numeric::odeint::integrate_adaptive(
                    boost::numeric::odeint::make_controlled<boost_stepper>( abs_error, rel_error ),
                    flow_equation, flow_state, 0.0, l_final, dl, boost::ref(book_keeper));
    }
    catch (LargeRODException& e) {
        std::cout << e.what() << std::endl;
    }

    // Checks whether symmetries are preserved
    {
        const auto& state = book_keeper.lowest_ROD_state;
        for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
            if (!NickelCUT::float_equal(state.dispersion[p], state.dispersion[-p])) {
                std::cerr << "Dispersion is not inversion symmetric: " << state.dispersion[p] << "  " << state.dispersion[-p] << std::endl;
                break;
            }
        }
    
        bool error_found = false;
        for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end() && !error_found; ++p) {
            for (momentum_iterator<L> q = momentum_iterator<L>::begin(); q != momentum_iterator<L>::end() && !error_found; ++q) {
                for (momentum_iterator<L> r = momentum_iterator<L>::begin(); r != momentum_iterator<L>::end() && !error_found; ++r) {
                    if(!NickelCUT::float_equal(state.interactions_differing_spin(p, q, r), state.interactions_differing_spin(-p, -q, -r))) {
                        std::cerr << "Interaction is not inversion symmetric " 
                            << state.interactions_differing_spin(p, q, r) 
                            << "  " << state.interactions_differing_spin(-p, -q, -r)
                            << std::endl;
                        error_found = true;
                    }
                }
            }
        }

        error_found = false;
        for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end() && !error_found; ++p) {
            for (momentum_iterator<L> q = momentum_iterator<L>::begin(); q != momentum_iterator<L>::end() && !error_found; ++q) {
                for (momentum_iterator<L> r = momentum_iterator<L>::begin(); r != momentum_iterator<L>::end() && !error_found; ++r) {
                    if(!NickelCUT::float_equal(state.interactions_differing_spin(p, q, r), state.interactions_differing_spin(p+r, q-r, -r))) {
                        std::cerr << "Interaction is not Hermitian " 
                            << state.interactions_differing_spin(p, q, r) 
                            << "  " << state.interactions_differing_spin(p+r, q-r, -r)
                            << std::endl;
                        error_found = true;
                    }
                }
            }
        }
    }

    const nlohmann::json j_metadata = model.generate_meta_data_json();
    nlohmann::json j_flow_data = book_keeper;
    j_flow_data.merge_patch(j_metadata);
    nlohmann::json j_full_flow_state = book_keeper.lowest_ROD_state;
    j_full_flow_state.merge_patch(j_metadata);

    mrock::utility::save_string(j_flow_data.dump(4), output_folder + data_file_names::FLOW_STEPS);
    mrock::utility::save_string(j_full_flow_state.dump(4), output_folder + data_file_names::FULL_FLOW_STATE);
    serialize_flow_state(book_keeper.lowest_ROD_state, output_folder, data_file_names::LOWEST_ROD_STATE);
    serialize_flow_state(flow_state, output_folder, data_file_names::FINAL_FLOW_STATE);

    book_keeper.print_final();
    return 0;
}