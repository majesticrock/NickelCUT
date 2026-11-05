#include "../sources/shared/Model.hpp"
#include "../sources/shared/occupation_numbers.hpp"

#include <mrock/utility/InputFileReader.hpp>

#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string>

namespace {
void verify_filling(const std::string& config_name) {
    const std::filesystem::path config_path =
        std::filesystem::path(NICKELCUT_SOURCE_DIR) / "params" / config_name;
    mrock::utility::InputFileReader input(config_path.string());
    const NickelCUT::Model model(input);

    if (!std::isfinite(model.E_F)) {
        throw std::runtime_error("The computed chemical potential is not finite.");
    }
    const double computed_filling = NickelCUT::compute_occupation_numbers(model);
    if (std::abs(computed_filling - model.filling) > 1e-10) {
        throw std::runtime_error("The computed occupations do not match the target filling.");
    }
    const std::string model_info = model.info();
    for (const std::string parameter : {"filling=", "E_F=", "temperature="}) {
        if (model_info.find(parameter) == std::string::npos) {
            throw std::runtime_error("The model summary omits " + parameter);
        }
    }
}
}  // namespace

int main() {
    verify_filling("half_filled.config");
    verify_filling("local.config");
    verify_filling("long_queue.config");
    verify_filling("mean_field_basic.config");
    verify_filling("med_queue.config");
    return 0;
}
