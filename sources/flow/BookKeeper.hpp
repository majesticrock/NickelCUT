#pragma once

#include "FlowContainer.hpp"

#include <nlohmann/json.hpp>

#include <deque>
#include <cstddef>
#include <chrono>

namespace NickelCUT::flow {

struct BookKeeper {    
    std::deque<double> l_times;
    std::deque<double> residual_offdiagonalities;
    std::deque<double> max_interactions;

    /////////////////////////////////////////////////////////

    BookKeeper(const FlowContainer& initial_flow_state, double _band_width, double _dl, std::chrono::minutes::rep _max_runtime_duration);

    void print_final(const FlowContainer& x, double l);

    void operator()(const FlowContainer& x, double l);

private:
    using clock = std::chrono::high_resolution_clock;

    const double dl;
    const double max_dl;
    const double min_interaction_difference;
    const double initial_band_width;
    const std::chrono::minutes::rep max_runtime_duration;

    const clock::time_point begin;
    clock::time_point last;

    std::size_t current_idx;
};

void to_json(nlohmann::json& j, const BookKeeper& book_keeper) noexcept;

} // namespace NickelCUT::flow
