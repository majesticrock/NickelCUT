#pragma once

#include <cstddef>
#include <chrono>

namespace NickelCUT::flow {

struct BookKeeper {
    double lowest_ROD;
    double l_of_lowest_ROD;
    std::size_t index_of_lowest_ROD;
    
    BookKeeper(double initial_ROD);

    bool process_step(double current_l, double ROD);

    void print_final() const;

private:
    using clock = std::chrono::high_resolution_clock;

    const clock::time_point begin;
    clock::time_point last;

    std::size_t current_idx;
};

} // namespace NickelCUT::flow
