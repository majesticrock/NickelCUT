#pragma once

#include <mrock/symbolic_operators/WickTermCollector.hpp>

#include <array>

namespace NickelCUT::commute::Hamiltonian {
using namespace mrock::symbolic_operators;

void export_as_flow_equation(const std::array<WickTermCollector, 3> flow_coeffs);

} // namespace NickelCUT::commute