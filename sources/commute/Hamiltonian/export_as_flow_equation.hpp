#pragma once

#include "../../experimental/WickOrderedCollector.hpp"

#include <array>

namespace NickelCUT::commute::Hamiltonian {
using namespace mrock::symbolic_operators;

void export_as_flow_equation(const std::array<experimental::WickOrderedCollector, 3> flow_coeffs);

} // namespace NickelCUT::commute