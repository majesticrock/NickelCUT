#include "../sources/flow/numerical_setup.hpp"

#include "../sources/mean_field/Model.hpp"
#include "../sources/L.hpp"

#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>
#include <mrock/utility/OutputConvenience.hpp>
#include <nlohmann/json.hpp>

#include <string>

#ifndef OUTPUT_DATA_DIR
#define OUTPUT_DATA_DIR "build/" //"../../data/nickel_cut/"
#endif

using namespace NickelCUT;

int main(int /*argc*/, char** /*argv*/) {
    const std::string load_binary_dir = std::string(OUTPUT_DATA_DIR) 
        + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/") // ensures that OUTPUT_DATA_DIR ends in "/"
        + "binaries/";
    mean_field::Model model(load_binary_dir, flow::U_0, flow::tprime, flow::mu_0, flow::T);

    const std::string output_dir = std::string(OUTPUT_DATA_DIR) 
        + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/") // ensures that OUTPUT_DATA_DIR ends in "/"
        + model.data_dir_name();
    
    auto solver = mrock::utility::Selfconsistency::make_iterative<double>(&model, &model.deltas, 1e-6);//, 200, 1e-6
    solver.compute();

    std::array<double, N> Delta_SC; 
    std::array<double, N> Delta_AFM;
    std::array<double, N> Delta_CDW;
    std::array<double, N> epsilon_up;
    std::array<double, N> epsilon_down;

    for (std::size_t i = 0U; i < Delta_SC.size(); ++i) {
        Delta_SC[i] = model.Delta_SC(i);
        Delta_AFM[i] = 0.5 * (model.Delta_DW_up(i) - model.Delta_DW_down(i));
        Delta_CDW[i] = 0.5 * (model.Delta_DW_up(i) + model.Delta_DW_down(i));
        epsilon_up[i] = model.dispersion_up(i);
        epsilon_down[i] = model.dispersion_down(i);
    }

    nlohmann::json jData = {
        { "time", mrock::utility::time_stamp() },
        { "Delta_SC", Delta_SC },
        { "Delta_AFM", Delta_AFM },
        { "Delta_CDW", Delta_CDW },
        { "epsilon_up", epsilon_up },
        { "epsilon_down", epsilon_down }
    };

    mrock::utility::save_string(jData.dump(4), output_dir + "mean_field_solution.json.gz");

    return 0;
}