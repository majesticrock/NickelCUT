#pragma once

#include <mrock/symbolic_operators/Term.hpp>

#include <sstream>
#include <vector>

namespace NickelCUT::commute {
using namespace mrock::symbolic_operators;

std::vector<Term> commutator_of_cut(std::ostringstream& oss);

} // namespace NickelCUT::commute