#include "restructure_terms.hpp"
#include "../../experimental/WickOrderedCollector.hpp"

#include <algorithm>

namespace NickelCUT::commute::PairAnnihilation {
using namespace mrock::symbolic_operators;

void restructure_bilinear_term(experimental::WickOrderedTerm& term)
{
    term.redistribute_momenta(term.wick_expression.back().momentum, 'q');
}

void restructure_quartic_term(experimental::WickOrderedTerm& term)
{
    term.redistribute_momenta(term.wick_expression.front().momentum, 'p');
    term.redistribute_momenta(term.wick_expression[2].momentum, 'q', {'p'});
 
    std::sort(term.sums.momenta.begin(), term.sums.momenta.end());
}

void advanced_clean_up(experimental::WickOrderedCollector& terms)
{
    for (auto& term : terms) {
        if (term.wick_expression.size() == 2U) {
            restructure_bilinear_term(term);
        }
        if (term.wick_expression.size() == 4U) {
            restructure_quartic_term(term);
        }
    }
}

} // namespace NickelCUT::commute