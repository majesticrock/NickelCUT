#pragma once

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
        << "\\end{align}\n. Here, the $\\ell$-dependence of the coefficients is implied."
        << "Note that, at $\\ell=0$, the interaction is independent of the momenta." << std::endl;
        
    std::vector<Term> cut_commutator = commutator(H, cut_generator);
    clean_up(cut_commutator);

    std::cout << "\nThe commutator reads\n\\begin{align}\n[H, \\eta] =" << cut_commutator
        << "\\end{align}" << std::endl;
    
    return cut_commutator;
}

int main(){
    std::vector<Term> commutator = commutator_of_cut();

    return 0;
}