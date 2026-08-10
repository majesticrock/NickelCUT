#pragma once

#include <mrock/symbolic_operators/Term.hpp>

#include <vector>

namespace NickelCUT::commute {
using namespace mrock::symbolic_operators;

std::vector<Term> commutator_of_cut();

} // namespace NickelCUT::commute