#include "Sparse.hpp"

#include "../ExtractionContainer.hpp"
#include "../../shared/FlowContainer.hpp"

namespace NickelCUT::flow::BookKeeper {
Sparse::Sparse(const FlowContainer& initial_flow_state,
                       FlowEquation& flow_equation,
                       double _band_width,
                       double _dl,
                       std::chrono::minutes::rep _max_runtime_duration)
    : base(initial_flow_state, flow_equation, _band_width, _dl, _max_runtime_duration) 
{}

void Sparse::print_final(const FlowContainer& x, double l) {
    basic_wrap_up(x, l);
}

void Sparse::operator()(const FlowContainer& x, double l) {
    basic_processing(x, l);
}

void to_json(nlohmann::json& j, const Sparse& book_keeper) noexcept {
    j = book_keeper.convert_to_json();
    j.update({{"extracted_channels", ExtractionContainer(book_keeper.last_good_state)}});
}

}  // namespace NickelCUT::flow::BookKeeper
