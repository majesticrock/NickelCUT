#include "../sources/L.hpp"
#include "../sources/flow/occupation_numbers.hpp"
#include "../sources/flow/momentum_iterator.hpp"
#include "../sources/flow/FlowContainer.hpp"
#include "../sources/flow/Model.hpp"
#include "../sources/flow/FlowEquation.hpp"
#include "../sources/helper_functions.hpp"
#include "../sources/flow/BookKeeper.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <mrock/utility/OutputConvenience.hpp>
#include <nlohmann/json.hpp>

#include <cmath>
#include <iostream>
#include <filesystem>

#ifndef OUTPUT_DATA_DIR
#define OUTPUT_DATA_DIR "build/" //"../../data/nickel_cut/"
#endif

using namespace NickelCUT;
using namespace NickelCUT::flow;
using namespace boost::numeric::odeint;

constexpr double U = -0.5;
constexpr double _ROD_0 = (U < 0. ? -1. : 1.) * U * L;

constexpr double abs_error = 1e-6;
constexpr double rel_error = 1e-6;

constexpr double l_final   = (1.5 / _ROD_0);
constexpr double target_dl = (1. / (5. * L * _ROD_0));
constexpr double dl = target_dl / 50.;

namespace boost { namespace numeric { namespace odeint {
template<>
struct vector_space_norm_inf< NickelCUT::flow::FlowContainer >
{
    typedef double result_type;
    double operator()( const NickelCUT::flow::FlowContainer &p ) const
    {
        return p.norm_inf();
    }
};
} } }
typedef runge_kutta_fehlberg78<FlowContainer, double, FlowContainer, double, vector_space_algebra> boost_stepper;

void serialize_flow_state(const FlowContainer& state, const std::string& output_dir) {
    const std::string file = output_dir + "lowest_ROD_state.bin";
    std::ofstream ofs(file, std::ios::binary);
    if (ofs.good()) {
        boost::archive::binary_oarchive oa(ofs);
        oa << state;
    }
    else {
        throw std::runtime_error("Outputstream for " + file + " is bad!");
    }
}

int main(int /*argc*/, char** /*argv*/) {
    const Model model{ U, 0.0, 0.0, -1.0 };
    FlowContainer flow_state(model);
    const double filling = compute_occupation_numbers(model, flow_state);
    std::cout << "\nConstructed initial states. The filling of the system is " << filling << std::endl;

    FlowEquation flow_equation;
    BookKeeper book_keeper(flow_state, target_dl);

    try {
        integrate_adaptive(make_controlled<boost_stepper>( abs_error, rel_error ),
                    flow_equation, flow_state, 0.0, l_final, dl, boost::ref(book_keeper));
    }
    catch (LargeRODException& e) {
        std::cout << e.what() << std::endl;
    }

    // Checks whether inversion symmetry is preserved (it should be!)
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

    const std::string output_folder = std::string(OUTPUT_DATA_DIR) + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/");
    std::filesystem::create_directories(output_folder);

    nlohmann::json j_meta_data = {
        { "time", mrock::utility::time_stamp() },
        { "L", L },
        { "filling", filling },
        { "U_0", model.U_0 },
        { "tprime", model.tprime },
        { "chemical_potential", model.chemical_potential },
        { "beta", model.beta }
    };
    nlohmann::json j_flow_data = book_keeper;
    j_flow_data.merge_patch(j_meta_data);

    mrock::utility::save_string(j_flow_data.dump(4), output_folder + "test.json.gz");
    serialize_flow_state(book_keeper.lowest_ROD_state, output_folder);

    book_keeper.print_final();
    return 0;
}