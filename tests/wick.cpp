#include "../sources/experimental/WickOrderedCollector.hpp"

#include <mrock/symbolic_operators/Term.hpp>

using namespace mrock::symbolic_operators;

int main() {
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

    const std::vector<WickOperatorTemplate> wick_templates({ 
        WickOperatorTemplate({Num_Comparison}, Momentum(), OperatorType::Number) 
    });

    //experimental::WickOrderedCollector kinetic_result = experimental::wick_decompose(H_kin, wick_templates);
    //std::cout << "Kinetic term decomposition:\n\\begin{align*}" 
    //    << kinetic_result 
    //    << "\\end{align*}" << std::endl;

    experimental::WickOrderedCollector interaction_result = experimental::wick_decompose(H_int, wick_templates);
    clean_wick_ordered_terms(interaction_result);
    
    std::cout << "Interaction term decomposition:\n\\begin{align*}" 
        << interaction_result 
        << "\\end{align*}" << std::endl;

    return 0;
}