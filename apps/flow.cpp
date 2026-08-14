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
#include <list>
#include <iostream>
#include <filesystem>
#include <chrono>

#ifndef OUTPUT_DATA_DIR
#define OUTPUT_DATA_DIR "build/" //"../../data/nickel_cut/"
#endif

using namespace NickelCUT;
using namespace NickelCUT::flow;
using namespace boost::numeric::odeint;

//constexpr double abs_error = 1e-10;
//constexpr double rel_error = 1e-6;

constexpr double l_final = 2.5;
constexpr double dl_output = 0.05;
constexpr double dl_use = dl_output / 10;


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
    const Model model{ 0.5, 0.0, 0.0, -1.0 };
    FlowContainer flow_state(model);
    compute_occupation_numbers(model, flow_state);

    FlowEquation flow_equation;

    std::list<double> l_times = {0.};
    std::list<double> residual_offdiagonalities = { flow_state.residual_offdiagonality() };
    std::list<FlowContainer> flow_states = {flow_state};

    double l = 0.;
    BookKeeper book_keeper(residual_offdiagonalities.back());

    std::list<FlowContainer>::const_iterator lowest_ROD_it = flow_states.end();
    --lowest_ROD_it;

    while(l < l_final) {
        //integrate_adaptive(make_controlled< runge_kutta_cash_karp54< FlowContainer > >( abs_error, rel_error ),
        //            flow_equation, flow_state, l, l + dl_output, dl_use);

        runge_kutta4<FlowContainer, double, FlowContainer, double, vector_space_algebra> stepper;
        integrate_const(stepper, flow_equation, flow_state, l , l + dl_output , dl_use );

        if (flow_state.norm_inf() > 1e3 || flow_state.contains_nan_or_inf()) break;
        l += dl_output;
        l_times.push_back(l);
        flow_states.push_back(flow_state);
        residual_offdiagonalities.push_back(flow_state.residual_offdiagonality());

        if(book_keeper.process_step(l, residual_offdiagonalities.back()))
        {
            lowest_ROD_it = flow_states.end();
            --lowest_ROD_it;
        };

        if (residual_offdiagonalities.back() > 50 * residual_offdiagonalities.front()) break;
    }

    nlohmann::json jData = {
        { "time", mrock::utility::time_stamp() },
        { "l_times", l_times },
        { "flow_states", flow_states },
        { "number_of_data_points", l_times.size() },
        { "L", L },
        { "residual_offdiagonalities", residual_offdiagonalities },
        { "lowest_ROD", book_keeper.lowest_ROD },
        { "index_of_lowest_ROD", book_keeper.index_of_lowest_ROD },
        { "l_of_lowest_ROD", book_keeper.l_of_lowest_ROD }
    };

    // Checks whether inversion symmetry is preserved (it should be!)
    for (const auto& state : flow_states) {
        for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
            if (!NickelCUT::float_equal(state.dispersion[p], state.dispersion[-p])) {
                std::cerr << "Dispersions: " << state.dispersion[p] << "  " << state.dispersion[-p] << std::endl;
            }
        }
    }

    const std::string output_folder = std::string(OUTPUT_DATA_DIR) + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/");
    std::filesystem::create_directories(output_folder);

    mrock::utility::save_string(jData.dump(4), output_folder + "test.json.gz");
    serialize_flow_state(*lowest_ROD_it, output_folder);

    book_keeper.print_final();
    return 0;
}