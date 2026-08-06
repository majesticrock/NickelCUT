#include "../sources/experimental/WickOrderedCollector.hpp"

#include <mrock/symbolic_operators/Commutation>

#include <vector>
#include <iostream>

using namespace mrock::symbolic_operators;

std::vector<Term> commutator_of_cut() {
    using namespace mrock::symbolic_operators;

    const Operator c_q_up     = Operator(Momentum('q'), Index::SpinUp,   false);
    const Operator c_q_down   = Operator(Momentum('q'), Index::SpinDown, false);
    const Operator c_q_sigma  = Operator(Momentum('q'), Index::Sigma,    false);

    const Term H_kin(1, 
        Coefficient::RealInversionSymmetric("\\varepsilon", MomentumList(Momentum('q'))),
        SumContainer{MomentumSum{'q'}, IndexSum{ Index::Sigma }},
        std::vector<Operator>({
            c_q_sigma.hermitian_conjugate(),
            c_q_sigma
        })
    );

    const Term H_int(1, 
        Coefficient("U", MomentumList({Momentum('q'), Momentum('p'), Momentum('r')})),
        SumContainer{MomentumSum{'q', 'p', 'r'}, IndexSum{}},
        std::vector<Operator>({
            c_q_up.hermitian_conjugate(),
            c_q_down.with_momentum('p').hermitian_conjugate(),
            c_q_down.with_momentum(Momentum("p-r")),
            c_q_up.with_momentum(Momentum("q+r"))
        })
    );

    const std::vector<Term> H = {H_kin, H_int};


    const Term cut_generator(1, 
        Coefficient("\\alpha", MomentumList({Momentum('q'), Momentum('p'), Momentum('r')})),
        SumContainer{MomentumSum{'q', 'p', 'r'}, IndexSum{}},
        std::vector<Operator>({
            c_q_up.hermitian_conjugate(),
            c_q_down.with_momentum('p').hermitian_conjugate(),
            c_q_down.with_momentum(Momentum("p-r")),
            c_q_up.with_momentum(Momentum("q+r"))
        })
    );

    std::cout << "We work with the Hamiltonian \n\\begin{align}\nH =" << H 
        << "\\end{align}\n, and the generator of the CUT\n\\begin{align}\n\\eta =" << cut_generator
        << "\\end{align}\n. Here, the $\\ell$-dependence of the coefficients is implied. "
        << "Note that, at $\\ell=0$, the interaction is independent of the momenta. " << std::endl;
        
    std::vector<Term> cut_commutator = commutator(H, cut_generator);
    clean_up(cut_commutator);

    std::cout << "\nThe commutator reads\n\\begin{align}\n[H, \\eta] =" << cut_commutator
        << "\\end{align}" << std::endl;
    
    return cut_commutator;
}

int main(){
    const std::vector<Term> commutator = commutator_of_cut();
    const std::vector<WickOperatorTemplate> wick_templates({ 
        WickOperatorTemplate({Num_Comparison}, Momentum(), OperatorType::Number) 
    });

    experimental::WickOrderedCollector normal_ordered_result = experimental::wick_decompose(commutator, wick_templates);
    experimental::clean_wick_ordered_terms(normal_ordered_result);

    std::erase_if(normal_ordered_result.terms, [](const experimental::WickOrderedTerm& term){
        return term.wick_expression.size() > 4U;
    });

    for (auto& term : normal_ordered_result) {
        const Momentum current_momentum = term.wick_expression.front().momentum;
        if (current_momentum != Momentum('q')) {
            std::size_t i=0;
            MomentumSymbol::name_type transformer = current_momentum[i].name;
            while (!term.sums.momenta.is_summed_over(transformer) && ++i < current_momentum.size()) {
                transformer = current_momentum[i].name;
            }
            term.transform_momentum_sum(transformer, Momentum(), );
        }
    }

    std::cout << "After normal ordering with respect to the Fermi sea, we omit any contribution with more than 4 operators. "
        << "The result reads\n\\begin{align*}\n\t"
        << normal_ordered_result << "\\end{align*}" << std::endl;

    return 0;
}