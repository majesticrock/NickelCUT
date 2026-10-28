#pragma once
/**
 * This file is legacy for dense output of the extracted channels during the flow and based on the ROD.
 * However, I noticed that the ROD is a poor indictor of when to stop the flow.
 * The new version simply saves the very last state and interrupts the flow, when some interaction matrix element
 * grows too large. Thereby, it forgoes saving a lot of unneeded data.
 *
 * This file is kept in case I do need dense output later down the line.
 */
#include "base.hpp"
#include "ExtractionContainer.hpp"

#include <nlohmann/json.hpp>

#include <chrono>
#include <deque>

namespace NickelCUT::flow::BookKeeper {
struct Dense : private base {
    std::deque<ExtractionContainer> extracted_channels;

    /////////////////////////////////////////////////////////

    Dense(const FlowContainer& initial_flow_state,
          FlowEquation& flow_equation,
          double _band_width,
          double _dl,
          std::chrono::minutes::rep _max_runtime_duration);

    void print_final(const FlowContainer& x, double l);

    void operator()(const FlowContainer& x, double l);

    inline const FlowContainer& get_last_good_state() const noexcept {
        return last_good_state;
    }

    friend void to_json(nlohmann::json& j, const Dense& book_keeper) noexcept;
};

}  // namespace NickelCUT::flow::BookKeeper
