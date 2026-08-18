#include "commute_and_print.hpp"

//#include "extract_flow_coefficients.hpp"
#include "restructure_terms.hpp"
//#include "export_as_flow_equation.hpp"

#include "../commutator_of_cut.hpp"
#include "../../experimental/WickOrderedCollector.hpp"

#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/ExpectationValues>

#include <vector>
#include <list>

namespace NickelCUT::commute::PairAnnihilation
{
using namespace mrock::symbolic_operators;
using namespace NickelCUT::commute;

void commute_and_print(std::ostringstream& oss)
{
    std::vector<Term> commutator = commutator_pair_annihilation(oss);
    const std::vector<WickOperatorTemplate> wick_templates = get_wick_templates();
    const std::vector<std::unique_ptr<WickSymmetry>> symmetries = get_symmetries();

    experimental::WickOrderedCollector normal_ordered_pair = experimental::wick_decompose(commutator, wick_templates);
    experimental::clean_wick_ordered_terms(normal_ordered_pair, symmetries);

    std::erase_if(normal_ordered_pair.terms, [](const experimental::WickOrderedTerm& term){
        return term.wick_expression.size() > 4U;
    });

    advanced_clean_up(normal_ordered_pair);

    oss << "After normal ordering with respect to the Fermi sea, as before, we omit any contribution with more than 4 operators. "
        << "The result reads\n\\begin{align*}\n\t"
        << normal_ordered_pair 
        << "\\end{align*}" << std::endl;

    oss << "The unique types of Wick-ordered expressions are\n\\begin{align*}\n";
    std::list<experimental::WickOrderedExpression> unique_pairs;
    for(const auto& term : normal_ordered_pair) {
        if (!exists_in(unique_pairs, term.wick_expression)) {
            unique_pairs.push_back(term.wick_expression);
            oss << "\t&" << term.wick_expression << "\\\\\n";
        }
    }
    oss << "\\end{align*}" << std::endl;
}
} // namespace NickelCUT::commute::PairAnnihilation
