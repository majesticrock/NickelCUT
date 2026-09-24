#include "../sources/flow/numerical_setup.hpp"

#include "../sources/mean_field/Model.hpp"
#include "../sources/mean_field/TCFinder.hpp"
#include "../sources/L.hpp"

#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>
#include <mrock/utility/OutputConvenience.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <deque>
#include <stdexcept>
#include <string>

#ifndef OUTPUT_DATA_DIR
#define OUTPUT_DATA_DIR "build/" //"../../data/nickel_cut/"
#endif

using namespace NickelCUT;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Invalid number of arguments: Use <path_to_executable> <configfile>" << std::endl;
        return -1;
    }
    mrock::utility::InputFileReader input(argv[1]);

    const std::string load_binary_dir = std::string(OUTPUT_DATA_DIR) 
        + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/") // ensures that OUTPUT_DATA_DIR ends in "/"
        + input.getString("output_dir") + "/"
        + "binaries/";
    mean_field::Model model(load_binary_dir, input);

    const std::string output_dir = std::string(OUTPUT_DATA_DIR) 
        + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/") // ensures that OUTPUT_DATA_DIR ends in "/"
        + input.getString("output_dir") + "/"
        + model.data_dir_name();
    
    const nlohmann::json j_metadata = model.generate_meta_data_json();

    auto solver = mrock::utility::Selfconsistency::make_broyden<double>(&model, &model.deltas, 200, 1e-6);
    solver.compute(true, 600);
    if (!model.deltas.converged)
    {
        std::cout << "Initial mean-field calculations did not converge. Retrying..." << std::endl;
        solver.compute(true, 600);

        if(!model.deltas.converged) {
            throw std::runtime_error("Initial mean-field calculations did not converge even after retry.");
        }
    }

    nlohmann::json jData = model.selfconsistency_to_json();
    std::cout << "Finished initial mean-field calculations. Found "
        << to_string(model.order_type(1e-10)) 
        << " order. Proceeding with T_c..." << std::endl;

    mean_field::TCFinder tc(model, 0.05);
    std::deque<mean_field::TransitionData> transition_data = tc.compute();

    jData.update({
        { "time", mrock::utility::time_stamp() } ,
        { "transition_data", transition_data }
    });
    jData.update(j_metadata);

    mrock::utility::save_string(jData.dump(4), output_dir + "mean_field_solution.json.gz");

    return 0;
}