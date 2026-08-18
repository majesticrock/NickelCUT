#pragma once

#include "../L.hpp"

#include <array>

namespace NickelCUT::flow
{
    struct Model;
    struct FlowContainer;

    extern std::array<double, N> occupation_numbers;

    // Fills the occupation_numbers array and returns their sum (the filling of the system)
    double compute_occupation_numbers(const Model& model, const FlowContainer& init);
} // namespace NickelCUT::flow
