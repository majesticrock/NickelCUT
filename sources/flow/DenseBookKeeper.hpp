#pragma once
/**
 * This file is legacy for dense output of the extracted channels during the flow and based on the ROD.
 * However, I noticed that the ROD is a poor indictor of when to stop the flow.
 * The new version simply saves the very last state and interrupts the flow, when some interaction matrix element
 * grows too large. Thereby, it forgoes saving a lot of unneeded data.
 * 
 * This file is kept in case I do need dense output later down the line.
 */
#include "FlowContainer.hpp"
#include "ExtractionContainer.hpp"

#include <nlohmann/json.hpp>

#include <deque>
#include <cstddef>
#include <chrono>

namespace NickelCUT::flow {

struct DenseBookKeeper {
    std::deque<double> l_times;
    std::deque<double> residual_offdiagonalities;
    std::deque<double> max_interactions;

    std::deque<ExtractionContainer> extracted_channels;

    FlowContainer last_good_state;

    /////////////////////////////////////////////////////////

    DenseBookKeeper(const FlowContainer& initial_flow_state, double _band_width, double _dl, std::chrono::minutes::rep _max_runtime_duration);

    void print_final(const FlowContainer& x, double l);

    void operator()(const FlowContainer& x, double l);

private:
    using clock = std::chrono::high_resolution_clock;

    const double dl;
    const double initial_band_width;
    const std::chrono::minutes::rep max_runtime_duration;

    const clock::time_point begin;
    clock::time_point last;

    std::size_t current_idx;
};

void to_json(nlohmann::json& j, const DenseBookKeeper& book_keeper) noexcept;

} // namespace NickelCUT::flow
