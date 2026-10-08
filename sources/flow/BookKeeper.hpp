#pragma once

#include "FlowContainer.hpp"

#include <nlohmann/json.hpp>

#include <deque>
#include <cstddef>
#include <chrono>

namespace NickelCUT::flow {

struct FlowEquation;

struct BookKeeper {    
    std::deque<double> l_times;
    std::deque<double> residual_offdiagonalities;
    std::deque<double> derivative_residual_offdiagonalities;
    std::deque<double> max_interactions;

    // Numerical errors acrue over time. This variable saves the last state that is still "good".
    // "good" means that all required symmetries are still preserved.
    FlowContainer last_good_state;

    /////////////////////////////////////////////////////////

    BookKeeper(const FlowContainer& initial_flow_state, FlowEquation& flow_equation, double _band_width, double _dl, std::chrono::minutes::rep _max_runtime_duration);

    void print_final(const FlowContainer& x, double l);

    void operator()(const FlowContainer& x, double l);

private:
    using clock = std::chrono::high_resolution_clock;

    double compute_derivative_rod(const FlowContainer& state, double l);

    FlowEquation& flow_equation;
    const double dl;
    const double initial_band_width;
    const std::chrono::minutes::rep max_runtime_duration;

    const clock::time_point begin;
    clock::time_point last;

    std::size_t current_idx;
};

void to_json(nlohmann::json& j, const BookKeeper& book_keeper) noexcept;

} // namespace NickelCUT::flow
