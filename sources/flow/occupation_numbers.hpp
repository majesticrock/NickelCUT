#pragma once

#include "../L.hpp"

#include <array>

namespace NickelCUT::flow
{
    struct Model;
    struct FlowContainer;

    extern std::array<double, N> occupation_numbers;

    void compute_occupation_numbers(const Model& model, const FlowContainer& init);
} // namespace NickelCUT::flow
