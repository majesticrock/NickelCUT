#pragma once

#include "../L.hpp"

#include <array>

namespace NickelCUT {
struct Model;

extern std::array<double, N> occupation_numbers;  ///< Cache of the occupation numbers <n>
extern std::array<double, N> hole_numbers;        ///< Cache of the occupation numbers of the holes <h> = 1 - <n>

// Fills the occupation_numbers array and returns their sum (the filling of the system)
double compute_occupation_numbers(const Model& model);
}  // namespace NickelCUT
