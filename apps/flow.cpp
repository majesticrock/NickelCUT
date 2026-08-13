#include "../sources/flow/L.hpp"
#include "../sources/flow/momentum_iterator.hpp"
#include "../sources/flow/FlowContainer.hpp"
#include "../sources/flow/Model.hpp"
#include "../sources/flow/FlowEquation.hpp"
#include "../sources/helper_functions.hpp"

#include <boost/numeric/odeint.hpp>
#include <mrock/utility/OutputConvenience.hpp>
#include <nlohmann/json.hpp>

#include <cmath>
#include <list>
#include <iostream>
#include <filesystem>

#ifndef OUTPUT_DATA_DIR
#define OUTPUT_DATA_DIR "../../data/nickel_cut/"
#endif

using namespace NickelCUT::flow;
using namespace boost::numeric::odeint;

//constexpr double abs_error = 1e-10;
//constexpr double rel_error = 1e-6;

constexpr double l_final = 0.8;
constexpr double dl_output = 0.05;
constexpr double dl_use = dl_output / 20;

int main() {
    Model model{0.5, 0.0, 0.0, 50};
    FlowContainer flow_state(model);

    FlowEquation flow_equation;

    std::list<double> l_times = {0.};
    std::list<FlowContainer> flow_states = {flow_state};

    //FlowContainer dHdl = FlowContainer();

    double l = 0.;
    while(l < l_final) {
        //integrate_adaptive(make_controlled< runge_kutta_cash_karp54< FlowContainer > >( abs_error, rel_error ),
        //            flow_equation, flow_state, l, l + dl_output, dl_use);

        std::cout << "l=" << l << std::endl;
        runge_kutta4<FlowContainer> stepper;
        integrate_const(stepper, flow_equation, flow_state, l , l + dl_output , dl_use );
        if (flow_state.norm_inf() > 1e3 || flow_state.contains_nan_or_inf()) break;

        l += dl_output;
        l_times.push_back(l);
        flow_states.push_back(flow_state);
    }

    nlohmann::json jData = {
        { "time", mrock::utility::time_stamp() },
        { "l_times", l_times },
        { "flow_states", flow_states },
        { "number_of_data_points", l_times.size() },
        { "L", L }
    };

    //for (const auto& state : flow_states) {
    //    for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
    //    for (momentum_iterator<L> q = momentum_iterator<L>::begin(); q != momentum_iterator<L>::end(); ++q) {
    //    for (momentum_iterator<L> r = momentum_iterator<L>::begin(); r != momentum_iterator<L>::end(); ++r) {
    //        if (!float_equal( state.interactions_same_spin(p, q, r),
    //                    state.interactions_same_spin(q, p, -r)
    //        )) {
    //            std::cerr << "Same spin: " << state.interactions_same_spin(p, q, r) << " vs. " 
    //                << state.interactions_same_spin(q, p, -r) << std::endl;
    //        }
//
    //        if (!float_equal( state.interactions_differing_spin(p, q, r),
    //                    state.interactions_differing_spin(q, p, -r)
    //        )) {
    //            std::cerr << "Differing spin: " << state.interactions_differing_spin(p, q, r) << " vs. " 
    //                << state.interactions_differing_spin(q, p, -r) << std::endl;
    //        }
    //    }}}
    //}

    for (const auto& state : flow_states) {
        for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
            if (!NickelCUT::float_equal(state.dispersion[p], state.dispersion[-p])) {
                std::cerr << "Dispersions: " << state.dispersion[p] << "  " << state.dispersion[-p] << std::endl;
            }
        }
    }

    //const std::string output_folder = std::string(OUTPUT_DATA_DIR) + "/";
    //std::filesystem::create_directories(output_folder);

    mrock::utility::save_string(jData.dump(4), "build/test.json.gz");
    return 0;
}