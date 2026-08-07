#pragma once

#include "experimental/WickOrderedCollector.hpp"

#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/WickSymmetry.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>

#include <cstddef>
#include <algorithm>
#include <array>
#include <vector>
#include <iostream>

namespace NickelCUT::commute {
using namespace mrock::symbolic_operators;

std::vector<Term> commutator_of_cut() {
    using namespace mrock::symbolic_operators;

    const Operator c_q_sigma        = Operator(Momentum('q'), Index::Sigma,      false);
    const Operator c_q_sigma_prime  = Operator(Momentum('q'), Index::SigmaPrime, false);

    const Term H_kin(1, 
        Coefficient::RealInversionSymmetric("\\varepsilon", MomentumList(Momentum('q'))),
        SumContainer{MomentumSum{'q'}, IndexSum{ Index::Sigma }},
        std::vector<Operator>({
            c_q_sigma.hermitian_conjugate(),
            c_q_sigma
        })
    );

    const Term H_int(1, 
        Coefficient("U", MomentumList({Momentum('q'), Momentum('p'), Momentum('r')}), 
            IndexWrapper{Index::Sigma, Index::SigmaPrime},
            false, true
        ),
        SumContainer{MomentumSum{'q', 'p', 'r'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
        std::vector<Operator>({
            c_q_sigma.hermitian_conjugate(),
            c_q_sigma_prime.with_momentum('p').hermitian_conjugate(),
            c_q_sigma_prime.with_momentum(Momentum("p-r")),
            c_q_sigma.with_momentum(Momentum("q+r"))
        })
    );

    const std::vector<Term> H = {H_kin, H_int};


    const Term cut_generator_half(1, 
        Coefficient("\\alpha", 
            MomentumList({Momentum('q'), Momentum('p'), Momentum('r')}), 
            IndexWrapper{Index::Sigma, Index::SigmaPrime},
            false, true
        ),
        SumContainer{MomentumSum{'q', 'p', 'r'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
        std::vector<Operator>({
            c_q_sigma.hermitian_conjugate(),
            c_q_sigma_prime.with_momentum('p').hermitian_conjugate(),
            c_q_sigma_prime.with_momentum(Momentum("p-r")),
            c_q_sigma.with_momentum(Momentum("q+r"))
        })
    );
    Term cut_generator_second_half = cut_generator_half.hermitian_conjugate();
    cut_generator_second_half.flip_sign();

    const std::vector<Term> cut_generator = {cut_generator_half, cut_generator_second_half};

    std::cout << "We work with the Hamiltonian \n\\begin{align}\nH =" << H 
        << "\\end{align}\n, and the generator of the CUT\n\\begin{align}\n\\eta =" << cut_generator
        << "\\end{align}\n. Here, the $\\ell$-dependence of the coefficients is implied. "
        << "Note that, at $\\ell=0$, the interaction is independent of the momenta. " << std::endl;
        
    std::vector<Term> cut_commutator = commutator(H, cut_generator);
    clean_up(cut_commutator);
    
    return cut_commutator;
}

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

    std::sort(term.operators.begin(), term.operators.end(), [](const WickOperator& l, const WickOperator& r) {
        if (l.type < r.type)
            return true;
        else if (l.type > r.type)
            return false;
        if (l.indizes.empty() || r.indizes.empty()) return false;
        return l.indizes[0] < r.indizes[0];
    });

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
    
}

void restructure_bilinear_term(experimental::WickOrderedTerm& term) {
    if(!term.is_bilinear()) {
        std::cerr << "Term must be bilinear!" << std::endl;
        abort();
    }

    // Rename indizes in the wick-ordered expression to one unified scheme
    term.redistribute_indizes(term.wick_expression[0].indizes[0], Index::Sigma);
    // Rename momenta in the wick-ordered expression to one unified scheme
    term.redistribute_momenta(term.wick_expression[0].momentum, 'q');

    restructure_expectation_values(term, Index::SigmaPrime, {Index::Sigma}, 'p', {'q'});

    std::sort(term.sums.momenta.begin(), term.sums.momenta.end());
    std::sort(term.sums.spins.begin(), term.sums.spins.end());
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
    term.redistribute_momenta(term.wick_expression[0].momentum, 'q');
    term.redistribute_momenta(term.wick_expression[1].momentum, 'p', {'q'});
    try {
        term.redistribute_momenta(term.wick_expression[2].momentum, 'r', {'q', 'p'});
        term.transform_momentum_sum('r', Momentum("?+p"), '?');
        term.rename_momenta('?', 'r');
        term.invert_momentum_sum('r');
    }
    catch (std::invalid_argument& e) {};

    restructure_expectation_values(term, Index::GeneralSpin_S, {Index::Sigma, Index::SigmaPrime}, 's', {'q', 'p', 'r'});

    std::sort(term.sums.momenta.begin(), term.sums.momenta.end());
    std::sort(term.sums.spins.begin(), term.sums.spins.end());
}

std::array<WickTermCollector, 3> extract_flow_coefficients(const experimental::WickOrderedCollector& terms) {
    std::array<WickTermCollector, 3> flow_coefficients;
    flow_coefficients[0].reserve(std::count_if(terms.begin(), terms.end(), 
        [](const experimental::WickOrderedTerm& term){
            return term.is_identity();
    }));
    flow_coefficients[1].reserve(std::count_if(terms.begin(), terms.end(), 
        [](const experimental::WickOrderedTerm& term){
            return term.is_bilinear();
    }));
    flow_coefficients[2].reserve(std::count_if(terms.begin(), terms.end(), 
        [](const experimental::WickOrderedTerm& term){
            return term.is_quartic();
    }));

    auto wick_term_from_normal_ordered = [](const experimental::WickOrderedTerm& term) {
        return WickTerm(
                term.multiplicity,
                term.coefficients,
                term.sums,
                term.delta_momenta,
                term.delta_indizes,
                term.operators
            );
    };

    const std::vector<Index> removeable_indizes_bilinear = {Index::Sigma};
    const std::vector<Index> removeable_indizes_quartic = {Index::Sigma, Index::SigmaPrime};

    const std::vector<MomentumSymbol::name_type> removeable_momenta_bilinear = {'q'};
    const std::vector<MomentumSymbol::name_type> removeable_momenta_quartic = {'q', 'p', 'r'};

    auto remove_sums = [](WickTerm& term, const std::vector<Index>& indizes, const std::vector<MomentumSymbol::name_type>& momenta) {
        auto& index_sums = term.sums.spins.summations;
        auto& momentum_sums = term.sums.momenta.summations;

        std::erase_if(index_sums, [&indizes](Index index) {
            return exists_in(indizes, index);
        });

        std::erase_if(momentum_sums, [&momenta](MomentumSymbol::name_type momentum) {
            return exists_in(momenta, momentum);
        });
    };

    for (const auto& term : terms) {
        if (term.is_identity()) {
            flow_coefficients[0].push_back(wick_term_from_normal_ordered(term));
        }
        else if(term.is_bilinear()) {
            flow_coefficients[1].push_back(wick_term_from_normal_ordered(term));
            remove_sums(flow_coefficients[1].back(), removeable_indizes_bilinear, removeable_momenta_bilinear);
        }
        else if (term.is_quartic()) {
            flow_coefficients[2].push_back(wick_term_from_normal_ordered(term));
            remove_sums(flow_coefficients[2].back(), removeable_indizes_quartic, removeable_momenta_quartic);
        }
    }

    return flow_coefficients;
} 

} // namespace NickelCUT::commute