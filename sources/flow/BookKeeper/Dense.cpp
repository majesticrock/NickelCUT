#include "Dense.hpp"

#include "../FlowEquation.hpp"

#include <deque>

namespace NickelCUT::flow::BookKeeper {
Dense::Dense(const FlowContainer& initial_flow_state,
                                 FlowEquation& flow_equation,
                                 double _band_width,
                                 double _dl,
                                 std::chrono::minutes::rep _max_runtime_duration)
    : base(initial_flow_state, flow_equation, _band_width, _dl, _max_runtime_duration),
    extracted_channels{ExtractionContainer(initial_flow_state)}
{}

void Dense::print_final(const FlowContainer& x, double l) {
    if (basic_wrap_up(x, l)) {
        extracted_channels.push_back(ExtractionContainer(x));
    }
}

void Dense::operator()(const FlowContainer& x, double l) {
    if (basic_processing(x, l)) {
        extracted_channels.push_back(ExtractionContainer(x));
    }
}

void to_json(nlohmann::json& j, const Dense& book_keeper) noexcept {
    j = book_keeper.convert_to_json();
    j.update({{"extracted_channels", book_keeper.extracted_channels}});
}

}  // namespace NickelCUT::flow::BookKeeper
