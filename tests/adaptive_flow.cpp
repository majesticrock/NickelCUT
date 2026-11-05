#include "../sources/flow/FlowEquation.hpp"
#include "../sources/shared/FlowContainer.hpp"
#include "../sources/shared/Model.hpp"
#include "../sources/shared/MomentumGrid.hpp"

#include <mrock/utility/InputFileReader.hpp>

#include <cmath>
#include <filesystem>
#include <string>

using namespace NickelCUT;

int main() {
    const std::filesystem::path config =
        std::filesystem::path(NICKELCUT_SOURCE_DIR) / "params" / "local.config";
    mrock::utility::InputFileReader input(config.string());
    const Model model(input);
    const auto fermion_count = MomentumGrid::fermion_points().size();
    const auto transfer_count = MomentumGrid::transfer_points().size();

    FlowContainer current(model);
    FlowContainer derivative;
    flow::FlowEquation{}(current, derivative, 0.0);

    if (fermion_count != MomentumGrid::fermion_points().size() ||
        transfer_count != MomentumGrid::transfer_points().size()) {
        return 1;
    }
    if (derivative.contains_nan_or_inf() || !std::isfinite(derivative.residual_offdiagonality())) {
        return 1;
    }
    if (!derivative.is_inversion_symmetric()) {
        return 1;
    }
    return 0;
}
