#pragma once

#include "../../shared/FlowContainer.hpp"

#include <nlohmann/json.hpp>

#include <chrono>
#include <cstddef>
#include <deque>

namespace NickelCUT::flow {
struct FlowEquation;
namespace BookKeeper {

class base {
protected:
    std::deque<double> l_times;
    std::deque<double> residual_offdiagonalities;
    std::deque<double> derivative_residual_offdiagonalities;
    std::deque<double> max_interactions;

    // Numerical errors acrue over time. This variable saves the last state that is still "good".
    // "good" means that all required symmetries are still preserved.
    FlowContainer last_good_state;

    /////////////////////////////////////////////////////////

    base(const FlowContainer& initial_flow_state,
               FlowEquation& flow_equation,
               double _band_width,
               double _dl,
               std::chrono::minutes::rep _max_runtime_duration);

    nlohmann::json convert_to_json() const noexcept;
    // Throws FlowExceptions to interrupt the flow if necessary and pushes the general data
    void basic_processing(const FlowContainer& x, const double l);
    // Returns true if the last step was accepted
    bool basic_wrap_up(const FlowContainer& x, const double l);

private:
    using clock = std::chrono::high_resolution_clock;

    double compute_derivative_rod(const FlowContainer& state, double l) const;

    FlowEquation& flow_equation;
    const double dl;
    const double initial_band_width;
    const std::chrono::minutes::rep max_runtime_duration;

    const clock::time_point begin;
    clock::time_point last;

    std::size_t current_idx;
};

} // namespace BookKeeper
} // namespace NickelCUT::flow
