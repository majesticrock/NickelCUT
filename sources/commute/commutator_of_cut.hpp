#pragma once

#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/ExpectationValues>

#include <sstream>
#include <memory>
#include <vector>

namespace NickelCUT::commute {
using namespace mrock::symbolic_operators;

TermCollector get_Hamiltonian();
TermCollector get_cut_generator();

std::vector<WickOperatorTemplate> get_wick_templates();
std::vector<std::unique_ptr<WickSymmetry>> get_symmetries();

TermCollector commutator_of_cut(std::ostringstream& oss);

TermCollector commutator_pair_annihilation(std::ostringstream& oss);

} // namespace NickelCUT::commute