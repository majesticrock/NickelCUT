#pragma once

#include "../../experimental/WickOrderedCollector.hpp"

namespace NickelCUT::commute::PairAnnihilation {
using namespace mrock::symbolic_operators;

void restructure_bilinear_term(experimental::WickOrderedTerm& term);

void restructure_quartic_term(experimental::WickOrderedTerm& term);

void advanced_clean_up(experimental::WickOrderedCollector& terms);

} // namespace NickelCUT::commute