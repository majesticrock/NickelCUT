#pragma once

#include "../L.hpp"

#include <array>

namespace NickelCUT::flow
{
    struct Model;

    extern std::array<double, N> occupation_numbers;

    // Fills the occupation_numbers array and returns their sum (the filling of the system)
    double compute_occupation_numbers(const Model& model);
} // namespace NickelCUT::flow
