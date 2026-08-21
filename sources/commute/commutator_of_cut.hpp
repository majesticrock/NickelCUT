#pragma once

#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/ExpectationValues>

#include <sstream>
#include <memory>
#include <vector>

namespace NickelCUT::commute {
using namespace mrock::symbolic_operators;

std::vector<Term> get_Hamiltonian();
std::vector<Term> get_cut_generator();

std::vector<WickOperatorTemplate> get_wick_templates();
std::vector<std::unique_ptr<WickSymmetry>> get_symmetries();

std::vector<Term> commutator_of_cut(std::ostringstream& oss);

std::vector<Term> commutator_pair_annihilation(std::ostringstream& oss);

} // namespace NickelCUT::commute