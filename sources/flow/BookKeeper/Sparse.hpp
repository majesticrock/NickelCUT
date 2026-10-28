#pragma once

#include "../FlowContainer.hpp"
#include "base.hpp"

#include <nlohmann/json.hpp>

#include <chrono>

namespace NickelCUT::flow::BookKeeper {
struct Sparse : private base {
    // No additional members

    Sparse(const FlowContainer& initial_flow_state,
               FlowEquation& flow_equation,
               double _band_width,
               double _dl,
               std::chrono::minutes::rep _max_runtime_duration);

    void print_final(const FlowContainer& x, double l);

    void operator()(const FlowContainer& x, double l);

    inline const FlowContainer& get_last_good_state() const noexcept {
        return last_good_state;
    }

    friend void to_json(nlohmann::json& j, const Sparse& book_keeper) noexcept;
};
}  // namespace NickelCUT::flow
