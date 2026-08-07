#include "../sources/experimental/WickOrderedCollector.hpp"
#include "../sources/commute_helper.hpp"

#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/WickSymmetry.hpp>

#include <memory>
#include <vector>
#include <list>
#include <iostream>

using namespace mrock::symbolic_operators;
using namespace NickelCUT::commute;

int main(){
    const std::vector<Term> commutator = commutator_of_cut();
    const std::vector<WickOperatorTemplate> wick_templates({ 
        WickOperatorTemplate({Num_Comparison}, Momentum(), OperatorType::Number) 
    });

    std::vector<std::unique_ptr<WickSymmetry>> symmetries;
    symmetries.push_back(std::make_unique<SpinSymmetry>());
    symmetries.push_back(std::make_unique<InversionSymmetry>());


    experimental::WickOrderedCollector normal_ordered_result = experimental::wick_decompose(commutator, wick_templates);
    experimental::clean_wick_ordered_terms(normal_ordered_result, symmetries);

    std::erase_if(normal_ordered_result.terms, [](const experimental::WickOrderedTerm& term){
        return term.wick_expression.size() > 4U;
    });

    for (auto& term : normal_ordered_result) {
        if (term.wick_expression.empty()) continue;

        // Rename indizes in the wick-ordered expression to one unified scheme
        term.redistribute_indizes(term.wick_expression[0].indizes[0], Index::Sigma);
        if (term.wick_expression.size() == 4U) {
            term.redistribute_indizes(term.wick_expression[1].indizes[0], Index::SigmaPrime, {Index::Sigma});

            // Assert that quartic terms are (sigma) (sigma') (sigma') (sigma)
            if (term.wick_expression[2].indizes[0] == Index::Sigma) {
                term.flip_sign();
                std::swap(term.wick_expression[2], term.wick_expression[3]);
            }
        }

        // Rename momenta in the wick-ordered expression to one unified scheme
        term.redistribute_momenta(term.wick_expression[0].momentum, 'q');

        if (term.wick_expression.size() == 4U) {
            term.redistribute_momenta(term.wick_expression[1].momentum, 'p', {'q'});
            
            try {
                term.redistribute_momenta(term.wick_expression[2].momentum, 'r', {'q', 'p'});
                term.transform_momentum_sum('r', Momentum("?+p"), '?');
                term.rename_momenta('?', 'r');
                term.invert_momentum_sum('r');
            }
            catch (std::invalid_argument& e) {};
        }
        std::sort(term.sums.momenta.begin(), term.sums.momenta.end());

        std::sort(term.sums.spins.begin(), term.sums.spins.end());
    }

    std::cout << "After normal ordering with respect to the Fermi sea, we omit any contribution with more than 4 operators. "
        << "The result reads\n\\begin{align*}\n\t"
        << normal_ordered_result << "\\end{align*}" << std::endl;

    std::cout << "The unique types of Wick-ordered expressions are\n\\begin{align*}\n";

    std::list<experimental::WickOrderedExpression> unique_wicks;
    for(const auto& term : normal_ordered_result) {
        if (!exists_in(unique_wicks, term.wick_expression)) {
            unique_wicks.push_back(term.wick_expression);
            std::cout << "\t&" << term.wick_expression << "\\\\\n";
        }
    }
    std::cout << "\\end{align*}" << std::endl;

    return 0;
}