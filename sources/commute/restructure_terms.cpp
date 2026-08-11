#include "restructure_terms.hpp"

#include <array>
#include <vector>

namespace NickelCUT::commute {
using namespace mrock::symbolic_operators;

void restructure_expectation_values(experimental::WickOrderedTerm& term, Index index_target,
        std::vector<Index> index_do_not_touch,
        MomentumSymbol::name_type momentum_target,
        std::vector<MomentumSymbol::name_type> momentum_do_not_touch) 
{
    for (auto& wick_operator : term.operators) {
        if (term.sums.spins.is_summed_over(wick_operator.indizes[0])) {
            while (exists_in(index_do_not_touch, index_target)) {
                ++index_target;
            }
            term.redistribute_indizes(wick_operator.indizes[0], index_target, index_do_not_touch);
            index_do_not_touch.push_back(index_target++);
        }
    }

    for (auto& wick_operator : term.operators) {
        while (exists_in(momentum_do_not_touch, momentum_target)) {
            ++(momentum_target._n);
        }
        try {
            term.redistribute_momenta(wick_operator.momentum, momentum_target, momentum_do_not_touch);
            momentum_do_not_touch.push_back(momentum_target);
            ++momentum_target._n;
        } catch (std::invalid_argument& e) { };
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
            } catch (std::invalid_argument& e) { };
        }
    }
    restructure_end_loop:

    term.sort();
}

void restructure_bilinear_term(experimental::WickOrderedTerm& term) {
    if(!term.is_bilinear()) {
        std::cerr << "Term must be bilinear!" << std::endl;
        abort();
    }

    // Rename indizes in the wick-ordered expression to one unified scheme
    term.redistribute_indizes(term.wick_expression[0].indizes[0], Index::Sigma);
    // Rename momenta in the wick-ordered expression to one unified scheme
    term.redistribute_momenta(term.wick_expression[0].momentum, 'p');

    restructure_expectation_values(term, Index::SigmaPrime, {Index::Sigma}, 'q', {'p'});
}

void restructure_quartic_term(experimental::WickOrderedTerm& term) {
    if(!term.is_quartic()) {
        std::cerr << "Term must be quartic!" << std::endl;
        abort();
    }

    // Rename indizes in the wick-ordered expression to one unified scheme
    term.redistribute_indizes(term.wick_expression[0].indizes[0], Index::Sigma);
    term.redistribute_indizes(term.wick_expression[1].indizes[0], Index::SigmaPrime, {Index::Sigma});
    // Assert that quartic terms are (sigma) (sigma') (sigma') (sigma)
    if (term.wick_expression[2].indizes[0] == Index::Sigma) {
        term.flip_sign();
        std::swap(term.wick_expression[2], term.wick_expression[3]);
    }

    // Rename momenta in the wick-ordered expression to one unified scheme
    term.redistribute_momenta(term.wick_expression[0].momentum, 'p');
    term.redistribute_momenta(term.wick_expression[1].momentum, 'q', {'p'});
    try {
        term.redistribute_momenta(term.wick_expression[2].momentum, 'r', {'p', 'q'});
        term.transform_momentum_sum('r', Momentum("?+q"), '?');
        term.rename_momenta('?', 'r');
        term.invert_momentum_sum('r');
    }
    catch (std::invalid_argument& e) {};

    restructure_expectation_values(term, Index::GeneralSpin_S, {Index::Sigma, Index::SigmaPrime}, 's', {'p', 'q', 'r'});
}

void advanced_clean_up(experimental::WickOrderedCollector& terms)
{
    for (auto& term : terms) {
        if (term.wick_expression.empty()) {
            restructure_expectation_values(term, Index::Sigma, {}, 'p', {});
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
        std::vector<MomentumSymbol::name_type>{'p'}, 
        std::vector<MomentumSymbol::name_type>{'p', 'q', 'r'} };

    for (auto& term : terms) {
        if (term.wick_expression.empty()) {
            restructure_coefficients(term, 'p', base_do_not_touch[0]);
        }
        else if (term.is_bilinear()) {
            restructure_coefficients(term, 'q', base_do_not_touch[1]);
        }
        else if (term.is_quartic()) {
            restructure_coefficients(term, 's', base_do_not_touch[2]);
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

    std::sort(terms.begin(), terms.end(), [](const experimental::WickOrderedTerm& l, const experimental::WickOrderedTerm& r){
        if (l.sums.momenta.size() < r.sums.momenta.size()) return true;
        if (l.sums.momenta.size() > r.sums.momenta.size()) return false;

        if (l.sums.spins.size() < r.sums.spins.size()) return true;
        if (l.sums.spins.size() > r.sums.spins.size()) return false;

        if (l.operators.size() < r.operators.size()) return true;
        if (l.operators.size() > r.operators.size()) return false;

        return false;
    });

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

        if (term.coefficients[0].momenta[2].is_used_at('r') >= 0 && term.coefficients[0].momenta[2].is_used_at('q') < 0) {
            term.rename_momenta('q', '_');
            term.rename_momenta('r', 'q');
            term.rename_momenta('_', 'r');
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

    //std::cerr << "# Terms: " << terms.size() << std::endl;
}

void improve_flow_coefficient_structure(WickTermCollector& terms)
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
            } catch (std::invalid_argument& e){};
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
            op.indizes[0] = Index::SpinUp;
        }

        std::sort(term.sums.momenta.begin(), term.sums.momenta.end());
    }

    terms.combine_duplicates();

    // Carry out the spin summations
    const std::size_t old_size = terms.size();
    terms.reserve(old_size + std::count_if(terms.begin(), terms.end(), [](const WickTerm& term) { return !term.sums.spins.empty();}));

    for (std::size_t i=0U; i<old_size; ++i){
        if (!terms[i].sums.spins.empty()) {
            assert(terms[i].sums.spins.size() == 1U);
            WickTerm copy = terms[i];
            
            for (std::size_t c=0U; c<terms[i].coefficients.size(); ++c) {
                Coefficient& coeff = terms[i].coefficients[c];
                if(coeff.indizes.size() < 2U) continue;

                terms[i].replace_each_index(terms[i].sums.spins[0], Index::SpinDown);
                copy.replace_each_index(terms[i].sums.spins[0], Index::SpinUp);
            }

            terms[i].sums.spins.clear();
        }
    }

    for (auto & term : terms) {
        for (auto& coeff : term.coefficients) {
            if (coeff.indizes.size() < 2U) continue;
            if (coeff.indizes[0] == coeff.indizes[1]) {
                coeff.indizes.clear();
                coeff.indizes.push_back(Index::Parallel);
            }
            else {
                coeff.indizes.clear();
                coeff.indizes.push_back(Index::AntiParallel);
            }
        }
    }

    std::sort(terms.begin(), terms.end(), [](const WickTerm& l, const WickTerm& r) {
        if (l.sums.momenta.size() < r.sums.momenta.size()) return true;
        if (l.sums.momenta.size() > r.sums.momenta.size()) return false;
        
        if (l.operators.size() < r.operators.size()) return true;
        if (l.operators.size() > r.operators.size()) return false;

        if (l.sums.spins.size() < r.sums.spins.size()) return true;
        if (l.sums.spins.size() > r.sums.spins.size()) return false;

        return false;
    });
}

}