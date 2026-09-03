#pragma once

#include "FlowContainer.hpp"
#include "DecouplingChannel.hpp"

#include <nlohmann/json.hpp>

#include <deque>
#include <cstddef>
#include <chrono>
#include <stdexcept>

namespace NickelCUT::flow {

class LargeRODException : public std::runtime_error {
public:
LargeRODException() : std::runtime_error("The ROD grew very large...") {};
};

struct ExtractionContainer {
    // The first element is for differing spins, the second for parallel spins
    std::pair<DecouplingChannel, DecouplingChannel> single_particle_energy;
    // The first element is for differing spins, the second for parallel spins
    std::pair<DecouplingChannel, DecouplingChannel> density_wave;
    // Only differing spins contribute
    DecouplingChannel superconductivity;

    FlowContainer::coeff_array dispersion;

    ExtractionContainer(const FlowContainer& x);
};

struct BookKeeper {
    double lowest_ROD;
    double l_of_lowest_ROD;
    std::size_t index_of_lowest_ROD;
    
    std::deque<double> l_times;
    std::deque<double> residual_offdiagonalities;

    std::deque<ExtractionContainer> extracted_channels;

    FlowContainer lowest_ROD_state;

    /////////////////////////////////////////////////////////

    BookKeeper(const FlowContainer& initial_flow_state, double _dl);

    // Returns true if the current ROD is the new lowest ROD
    bool process_step(double current_l, double ROD);

    void print_final() const;

    void operator()(const FlowContainer &x, double l);

private:
    using clock = std::chrono::high_resolution_clock;

    const double dl;
    const double max_dl;
    const double min_ROD_difference;

    const clock::time_point begin;
    clock::time_point last;

    std::size_t current_idx;
};

void to_json(nlohmann::json& j, const ExtractionContainer& extracted_channels) noexcept;

void to_json(nlohmann::json& j, const BookKeeper& book_keeper) noexcept;

} // namespace NickelCUT::flow
