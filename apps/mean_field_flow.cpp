#include "../sources/flow/ExtractionContainer.hpp"
#include "../sources/flow/Model.hpp"
#include "../sources/flow/data_file_names.hpp"
#include "../sources/mean_field/Model.hpp"
#include "../sources/mean_field/OrderType.hpp"
#include "../sources/L.hpp"

#include <mrock/utility/InputFileReader.hpp>
#include <mrock/utility/OutputConvenience.hpp>
#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>
#include <nlohmann/json.hpp>

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace NickelCUT;

#ifndef OUTPUT_DATA_DIR
#define OUTPUT_DATA_DIR "build/"
#endif

namespace
{
nlohmann::json load_flow_json(const std::filesystem::path& path)
{
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open flow data file: " + path.string());
    }

    boost::iostreams::filtering_istream decompressed;
    decompressed.push(boost::iostreams::gzip_decompressor());
    decompressed.push(file);
    return nlohmann::json::parse(decompressed);
}

std::string output_filename(const std::filesystem::path& flow_file)
{
    constexpr std::string_view flow_prefix = "flow.json.gz";
    const std::string filename = flow_file.filename().string();
    if (!filename.starts_with(flow_prefix)) {
        throw std::invalid_argument("Expected a flow.json.gz[segment] input file");
    }

    const std::string suffix = filename.substr(flow_prefix.size());
    if (!std::all_of(suffix.begin(), suffix.end(), [](unsigned char character) {
            return std::isdigit(character);
        })) {
        throw std::invalid_argument("Flow segment suffix must contain only digits");
    }

    return "mean_field_across_flow" + (suffix.empty() ? "" : "_" + suffix) + ".json.gz";
}
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <configfile>" << std::endl;
        return 1;
    }

    try {
        mrock::utility::InputFileReader input(argv[1]);
        const flow::Model model_parameters(input);
        const std::filesystem::path output_root = std::filesystem::path(OUTPUT_DATA_DIR)
            / input.getString("output_dir");
        const std::filesystem::path flow_file = output_root
            / model_parameters.data_dir_name()
            / flow::data_file_names::FLOW_STEPS;

        if (!std::filesystem::is_regular_file(flow_file)) {
            throw std::runtime_error("Flow data file not found: " + flow_file.string());
        }

        const nlohmann::json flow_data = load_flow_json(flow_file);
        if (flow_data.at("L").get<int>() != L) {
            throw std::invalid_argument("Flow data lattice size does not match this executable");
        }

        const auto& l_times_json = flow_data.at("l_times");
        const auto& rods_json = flow_data.at("residual_offdiagonalities");
        const auto& extracted_json = flow_data.at("extracted_channels");
        if (!l_times_json.is_array() || !rods_json.is_array() || !extracted_json.is_array()
            || l_times_json.empty() || l_times_json.size() != rods_json.size()
            || l_times_json.size() != extracted_json.size()) {
            throw std::invalid_argument("Flow data step arrays are missing or have inconsistent sizes");
        }

        const std::vector<double> l_times = l_times_json.get<std::vector<double>>();
        const std::vector<double> rods = rods_json.get<std::vector<double>>();
        const double flow_temperature = flow_data.at("T").get<double>();
        mean_field::Model model(
            extracted_json.at(0).get<flow::ExtractionContainer>(),
            input
        );

        model.temperature = 0.0;
        model.beta = -1.0;

        nlohmann::json output = {
            { "time", mrock::utility::time_stamp() },
            { "source_flow_file", flow_file.string() },
            { "L", flow_data.at("L") },
            { "flow_temperature", flow_temperature },
            { "mean_field_temperature", 0.0 },
            { "U_0", flow_data.at("U_0") },
            { "tprime", flow_data.at("tprime") },
            { "E_F", flow_data.at("E_F") },
            { "target_filling", model.target_filling },
            { "l_times", l_times_json },
            { "residual_offdiagonalities", rods_json },
            { "solutions", nlohmann::json::array() }
        };

        for (std::size_t index = 0; index < l_times.size(); ++index) {
            model.extracted_channels = extracted_json.at(index).get<flow::ExtractionContainer>();
            model.reset_self_consistency_values();
            model.deltas.converged = false;

            auto solver = mrock::utility::Selfconsistency::make_broyden<double>(
                &model, &model.deltas, 200, 1e-6
            );
            solver.compute(true, 600);
            if (!model.deltas.converged) {
                std::cout << "Mean-field calculation did not converge at l=" << l_times[index]
                    << ". Retrying..." << std::endl;
                solver.compute(true, 600);
                if (!model.deltas.converged) {
                    throw std::runtime_error("Mean-field calculation did not converge at l="
                        + std::to_string(l_times[index]));
                }
            }

            nlohmann::json solution = model.selfconsistency_to_json();
            solution.update({
                { "l", l_times[index] },
                { "residual_offdiagonality", rods[index] },
                { "chemical_potential", model.chemical_potential },
                { "order_type", to_string(model.order_type(1e-10)) }
            });
            output["solutions"].push_back(std::move(solution));
            std::cout << "Finished zero-temperature mean-field calculation at l="
                << l_times[index] << " (" << index + 1 << "/" << l_times.size() << ")"
                << std::endl;
        }

        const std::filesystem::path output_path = flow_file.parent_path() / output_filename(flow_file);
        mrock::utility::save_string(output.dump(4), output_path.string());
        std::cout << "Saved results to " << output_path << std::endl;
    }
    catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }

    return 0;
}