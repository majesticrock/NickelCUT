#include "restructure_terms.hpp"
#include "../../experimental/WickOrderedCollector.hpp"

#include <mrock/symbolic_operators/detail/OperatorOrder.hpp>

#include <array>
#include <vector>
#include <list>
#include <string>
#include <utility>

namespace NickelCUT::commute::Hamiltonian {
using namespace mrock::symbolic_operators;

void restructure_expectation_values(experimental::WickOrderedTerm& term,
        MomentumSymbol::name_type momentum_target,
        std::vector<MomentumSymbol::name_type> momentum_do_not_touch) 
{
    for (auto& wick_operator : term.operators) {
        wick_operator.indices.front() = Index::SpinUp; // Use the spin symmetry
    }

    for (auto& wick_operator : term.operators) {
        while (exists_in(momentum_do_not_touch, momentum_target)) {
            ++(momentum_target._n);
        }
        try {
            term.redistribute_momenta(wick_operator.momentum, momentum_target, momentum_do_not_touch);
            momentum_do_not_touch.push_back(momentum_target);
            ++momentum_target._n;
        } catch (redistribution_error& e) { };
    }

    term.sort();
}

void restructure_coefficients(experimental::WickOrderedTerm& term,
        MomentumSymbol::name_type momentum_target,
        std::vector<MomentumSymbol::name_type> momentum_do_not_touch)
{
    for (auto& coeff : term.coefficients) {
        for (auto& coeff_momentum : coeff.momenta) {
            try {
                term.redistribute_momenta(coeff_momentum, momentum_target, momentum_do_not_touch);
                momentum_do_not_touch.push_back(momentum_target);
                ++momentum_target._n;
                if (momentum_do_not_touch.size() >= term.sums.momenta.size()) goto restructure_end_loop;
            } catch (redistribution_error& e) { };
        }
    }
    restructure_end_loop:

    for (auto& coeff : term.coefficients) {
        coeff.use_custom_symmetry();
    }

    term.sort();
}

void restructure_bilinear_term(experimental::WickOrderedTerm& term) {
    if(!term.is_bilinear()) {
        std::cerr << "Term must be bilinear!" << std::endl;
        abort();
    }

    // Rename indices in the wick-ordered expression to one unified scheme
    term.redistribute_indices(term.wick_expression[0].indices[0], Index::Sigma);

    assert(term.wick_expression[0].indices[0] == Index::Sigma);

    restructure_expectation_values(term, 'P', {'K'});
}

void restructure_quartic_term(experimental::WickOrderedTerm& term) {
    if(!term.is_quartic()) {
        std::cerr << "Term must be quartic!" << std::endl;
        abort();
    }

    // Rename indices in the wick-ordered expression to one unified scheme
    term.redistribute_indices(term.wick_expression[0].indices[0], Index::Sigma);
    term.redistribute_indices(term.wick_expression[1].indices[0], Index::SigmaPrime, {Index::Sigma});
    
    assert(term.wick_expression[0].indices[0] == Index::Sigma);
    assert(term.wick_expression[1].indices[0] == Index::SigmaPrime);
    assert(term.wick_expression[2].indices[0] == Index::SigmaPrime);
    assert(term.wick_expression[3].indices[0] == Index::Sigma);

    restructure_expectation_values(term, 'R', {'K', 'P', 'Q'});
}

void advanced_clean_up(experimental::WickOrderedCollector& terms)
{
    {
        auto it = std::remove_if(terms.begin(), terms.end(), [](const experimental::WickOrderedTerm& term){
            for (const auto& coeff : term.coefficients) {
                if (coeff.name != "\\alpha") continue;
                // The alpha-coefficient is proportional to sign(eps(k) + eps(p) - eps(p-q) - eps(k+q))
                // Now, if q=0 (the last Momentum in momenta), then this reads sign(0), which is just 0
                // Hence, the term does not contribute
                if (coeff.momenta.back() == Momentum()) return true;
            }
            return false;
        });
        terms.erase(it, terms.end());
    }

    for (auto& term : terms) {
        if (term.wick_expression.empty()) {
            restructure_expectation_values(term, 'K', {});
        }
        else if (term.is_bilinear()) {
            restructure_bilinear_term(term);
        }
        else if (term.is_quartic()) {
            restructure_quartic_term(term);
        }
        for (auto& coeff : term.coefficients) {
            if (coeff.momenta.size() == 3U){
                if (coeff.momenta[0] > coeff.momenta[1]) {
                    coeff.use_symmetric_interaction_exchange();
                }
                if (coeff.momenta[0].first_momentum_is_negative()) {
                    coeff.use_symmetric_interaction_inversion();
                }
                coeff.use_custom_symmetry();
            }
            else if (coeff.momenta.size() == 1U) {
                if (coeff.momenta[0].first_momentum_is_negative()) {
                    coeff.use_inversion_symmetry();
                }
            }
        }
        for (auto& op : term.operators) {
            if (op.momentum.first_momentum_is_negative()) {
                // Inversion symmetry
                op.momentum.flip_momentum();
            }
        }
    }
    terms.combine_duplicates();

    const std::array<std::vector<MomentumSymbol::name_type>, 3> base_do_not_touch{ std::vector<MomentumSymbol::name_type>{}, 
        std::vector<MomentumSymbol::name_type>{'K'}, 
        std::vector<MomentumSymbol::name_type>{'K', 'P', 'Q'} };

    for (auto& term : terms) {
        if (term.wick_expression.empty()) {
            restructure_coefficients(term, 'K', base_do_not_touch[0]);
        }
        else if (term.is_bilinear()) {
            restructure_coefficients(term, 'P', base_do_not_touch[1]);
        }
        else if (term.is_quartic()) {
            restructure_coefficients(term, 'R', base_do_not_touch[2]);
        }

        for (auto& coeff : term.coefficients) {
            if (coeff.momenta.size() == 3U) {
                if (!coeff.momenta[0].empty() && !coeff.momenta[1].empty()) {
                    if (coeff.momenta[0].front() > coeff.momenta[1].front()) {
                        coeff.use_symmetric_interaction_exchange();
                    }
                } 
                if (coeff.momenta[0].first_momentum_is_negative()) {
                    coeff.use_symmetric_interaction_inversion();
                }
            }
        }

        for (auto& op : term.operators) {
            if (op.momentum.first_momentum_is_negative()) {
                // Inversion symmetry
                op.momentum.flip_momentum();
            }
        }
    }
    terms.combine_duplicates();

    for (auto& term : terms) {
        std::vector<MomentumSymbol::name_type> do_not_touch = base_do_not_touch[term.wick_expression.size() / 2];
        for (auto& summed : term.sums.momenta) {
            if (exists_in(do_not_touch, summed)) continue;
            if (term.coefficients[0].momenta.factor_of_first_occurance(summed) < 0) {
                term.invert_momentum_sum(summed);
                do_not_touch.push_back(summed);
            }
        }

        for (auto& coeff : term.coefficients) {
            if (coeff.momenta[0] > coeff.momenta[1]) {
                coeff.use_symmetric_interaction_exchange();
            }
            if (coeff.momenta[0].first_momentum_is_negative()) {
                coeff.use_symmetric_interaction_inversion();
            }
        }

        for (auto& op : term.operators) {
            if (op.momentum.first_momentum_is_negative()) {
                // Inversion symmetry
                op.momentum.flip_momentum();
            }
        }
    }
    terms.combine_duplicates();

    for (auto& term : terms) {
        // Special treatment for bilinears
        if (!term.is_bilinear()) continue;
        if (term.coefficients[0].momenta.size() != 3U) continue;

        if (term.coefficients[0].momenta[2].is_used_at('Q') >= 0 && term.coefficients[0].momenta[2].is_used_at('P') < 0) {
            term.rename_momenta('P', PLACEHOLDER_SYMBOL);
            term.rename_momenta('Q', 'P');
            term.rename_momenta(PLACEHOLDER_SYMBOL, 'Q');
        }

        for (auto& coeff : term.coefficients) {
            if (coeff.momenta[0] > coeff.momenta[1]) {
                coeff.use_symmetric_interaction_exchange();
            }
            if (coeff.momenta[0].first_momentum_is_negative()) {
                coeff.use_symmetric_interaction_inversion();
            }
        }
        
        term.sort();
    }
    terms.combine_duplicates();

    // Bonus help for bilinear Hartree-like terms
    for (auto& term : terms) {
        if (!term.is_bilinear() || term.coefficients[0].momenta.size() != 3U) continue;
        if (term.coefficients[0].momenta.back() != Momentum()) continue;
        if (term.coefficients[1].momenta[1] == Momentum('P')) {
            term.swap_momenta('P', 'Q');
            term.coefficients[0].use_symmetric_interaction_exchange();
            if (term.operators.size() == 2U) 
                std::swap(term.operators[0], term.operators[1]);
        }
    }
    // Bonus help for quartic Hartree-like terms 
    for (auto& term : terms) {
        if (!term.is_quartic() || term.coefficients[0].momenta.size() != 3U) continue;
        if (term.coefficients[0].momenta.back() != Momentum()) continue;
        if (term.coefficients[0].momenta[1] == Momentum("K+Q")) {
            term.swap_momenta('K', 'P');
            term.invert_momentum_sum('Q');

            term.swap_indices(Index::Sigma, Index::SigmaPrime);
            sort_operators_by_indices(term.wick_expression.operators);
            for (auto& coeff : term.coefficients) {
                coeff.use_custom_symmetry();
            }

            term.coefficients[1].use_symmetric_interaction_exchange();
        }
    }

    terms.sort();
    terms.combine_duplicates();

    std::cout << "# Terms: " << terms.size() << std::endl;
}

void improve_flow_coefficient_structure(experimental::WickOrderedCollector& terms)
{
    auto is_transformable = [](const WickOperator& op, const MomentumSum& sums) {
        for (const auto& symbol : op.momentum) {
            if (sums.is_summed_over(symbol.name)) return true;
        }
        return false;
    };

    for (auto& term : terms) {
        if (term.sums.momenta.empty()) continue;
        std::vector<MomentumSymbol::name_type> used;
        MomentumSymbol::name_type target = term.sums.momenta.front();
        for (auto& op : term.operators) {
            if (!is_transformable(op, term.sums.momenta)) continue;
            try {
                term.redistribute_momenta(op.momentum, target, used);
                used.push_back(target);
                ++target._n;
            } catch (redistribution_error& e){};
        }

        for (auto& coeff : term.coefficients) {
            if (coeff.momenta.size() == 3U) {
                if (coeff.momenta[0] > coeff.momenta[1]) {
                    coeff.use_symmetric_interaction_exchange();
                }
                if (coeff.momenta.back().first_momentum_is_negative()) {
                    coeff.use_symmetric_interaction_inversion();
                }
            }
            
        }

        for (auto& op : term.operators) {
            // Fermi sea is independent of the spin orientation
            op.indices[0] = Index::SpinUp;
        }
        term.sort();
    }

    terms.combine_duplicates();

    // Carry out the spin summations
    const std::size_t old_size = terms.size();
    terms.reserve(old_size + std::count_if(terms.begin(), terms.end(), [](const experimental::WickOrderedTerm& term) { 
        return !term.sums.spins.empty();
    }));

    // This runs the SigmaPrime sum for bilinears and the GeneralSpin_S sum for quartics
    for (std::size_t i=0U; i<old_size; ++i){
        if (!terms[i].sums.spins.empty()) {
            assert(terms[i].sums.spins.size() == 1U);
            experimental::WickOrderedTerm copy = terms[i];
            
            terms[i].replace_each_index(terms[i].sums.spins[0], Index::SpinDown);
            copy.replace_each_index(copy.sums.spins[0], Index::SpinUp);

            terms[i].sums.spins.clear();
            copy.sums.spins.clear();

            terms.push_back(copy);
        }
    }

    for (auto & term : terms) {
        // We simply set Sigma=SpinUp; the entire problem is invariant under SpinUp <-> SpinDown
        // Hence, the spin down states will transform exactly the same as the spin up ones.
        // Thus, replacing Sigma by either direction in the coefficients does not change anything
        // It is important to remember though that the Sigma summation for :c_(k,sigma)^dagger c_(k,sigma):
        // does remain! We only abuse that U_(sigma sigma) = U_(down down)
        for (auto& coeff : term.coefficients) {
            if (coeff.indices.size() < 2U) continue;
            coeff.indices.replace_index(Index::Sigma, Index::SpinUp);
            coeff.indices.replace_index(Index::AntiSigma, Index::SpinDown);
            
            assert(coeff.indices[0] == Index::SpinDown || coeff.indices[0] == Index::SpinUp);
            assert(coeff.indices[1] == Index::SpinDown || coeff.indices[1] == Index::SpinUp);

            if (coeff.indices[0] == coeff.indices[1]) {
                coeff.indices.clear();
                coeff.indices.push_back(Index::Parallel);
            }
            else {
                coeff.indices.clear();
                coeff.indices.push_back(Index::AntiParallel);
            }
        }
    }

    terms.combine_duplicates();
    terms.sort();
}

}