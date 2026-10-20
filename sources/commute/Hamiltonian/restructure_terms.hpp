#pragma once

#include <mrock/symbolic_operators/IndexWrapper.hpp>
#include <mrock/symbolic_operators/MomentumSymbol.hpp>
#include <mrock/symbolic_operators/experimental/WickOrderedCollector.hpp>
#include <mrock/symbolic_operators/experimental/WickOrderedTerm.hpp>

#include <vector>

namespace NickelCUT::commute::Hamiltonian {
using namespace mrock::symbolic_operators;

void restructure_expectation_values(experimental::WickOrderedTerm& term,
                                    MomentumSymbol::name_type momentum_target,
                                    std::vector<MomentumSymbol::name_type> momentum_do_not_touch);

void restructure_coefficients(experimental::WickOrderedTerm& term,
                              MomentumSymbol::name_type momentum_target,
                              std::vector<MomentumSymbol::name_type> momentum_do_not_touch);

void restructure_bilinear_term(experimental::WickOrderedTerm& term);

void restructure_quartic_term(experimental::WickOrderedTerm& term);

void advanced_clean_up(experimental::WickOrderedCollector& terms);

void improve_flow_coefficient_structure(experimental::WickOrderedCollector& terms);

}  // namespace NickelCUT::commute::Hamiltonian