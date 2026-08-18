#pragma once

#include "../../experimental/WickOrderedCollector.hpp"

#include <mrock/symbolic_operators/WickTerm.hpp>

#include <array>

namespace NickelCUT::commute::Hamiltonian {
using namespace mrock::symbolic_operators;

std::array<WickTermCollector, 3> extract_flow_coefficients(const experimental::WickOrderedCollector& terms);

} // namespace NickelCUT::commute