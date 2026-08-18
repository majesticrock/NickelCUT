#include "../sources/mean_field/Model.hpp"
#include "../sources/L.hpp"

#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>
#include <mrock/utility/OutputConvenience.hpp>
#include <nlohmann/json.hpp>

#include <string>

#ifndef OUTPUT_DATA_DIR
#define OUTPUT_DATA_DIR "build/" //"../../data/nickel_cut/"
#endif

int main(int /*argc*/, char** /*argv*/) {
    const std::string output_folder = std::string(OUTPUT_DATA_DIR) + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/");

    NickelCUT::mean_field::Model model(output_folder);

    auto solver = mrock::utility::Selfconsistency::make_broyden<double>(&model, &model.deltas, 200, 1e-6);
    solver.compute();

    std::array<double, NickelCUT::N> Delta_SC; 
    std::array<double, NickelCUT::N> Delta_AFM;
    std::array<double, NickelCUT::N> Delta_CDW;
    std::array<double, NickelCUT::N> epsilon_up;
    std::array<double, NickelCUT::N> epsilon_down;

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

    mrock::utility::save_string(jData.dump(4), output_folder + "mean_field_solution.json.gz");

    return 0;
}