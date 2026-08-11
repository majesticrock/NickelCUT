#include "commutator_of_cut.hpp"

#include "../experimental/WickOrderedCollector.hpp"
#include "restructure_terms.hpp"

#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/ExpectationValues>

namespace NickelCUT::commute {
using namespace mrock::symbolic_operators;

std::vector<Term> commutator_of_cut(std::ostringstream& oss) {
    using namespace mrock::symbolic_operators;

    const Operator c_p_sigma        = Operator(Momentum('p'), Index::Sigma,      false);
    const Operator c_p_sigma_prime  = Operator(Momentum('p'), Index::SigmaPrime, false);

    // It only matters if the indizes of the coefficients are identical or differ;
    // e.g. U_(sigma, sigma') = U_(simga', sigma), and U_(down, down) = U_(up up)
    const std::function<void(Coefficient&)> spin_symmetry = [](Coefficient& coeff) {
        assert(coeff.indizes.size() == 2U);
        if (coeff.indizes[0] > coeff.indizes[1]) {
            std::swap(coeff.indizes[0], coeff.indizes[1]);
        }
    };

    const Term H_kin(1, 
        Coefficient::RealInversionSymmetric("\\tilde{\\varepsilon}", MomentumList(Momentum('p'))),
        SumContainer{MomentumSum{'p'}, IndexSum{ Index::Sigma }},
        std::vector<Operator>({
            c_p_sigma.hermitian_conjugate(),
            c_p_sigma
        })
    );

    const Term H_int(1, 
        Coefficient::RealInteraction("U", MomentumList({Momentum('p'), Momentum('q'), Momentum('r')}), 
            IndexWrapper{Index::Sigma, Index::SigmaPrime},
            spin_symmetry
        ),
        SumContainer{MomentumSum{'p', 'q', 'r'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
        std::vector<Operator>({
            c_p_sigma.hermitian_conjugate(),
            c_p_sigma_prime.with_momentum('q').hermitian_conjugate(),
            c_p_sigma_prime.with_momentum(Momentum("q-r")),
            c_p_sigma.with_momentum(Momentum("p+r"))
        })
    );

    const std::vector<Term> H = {H_kin, H_int};

    // --------------------------------------------------------- //

    const std::vector<WickOperatorTemplate> wick_templates({ 
        WickOperatorTemplate({Num_Comparison}, Momentum(), OperatorType::Number) 
    });

    std::vector<std::unique_ptr<WickSymmetry>> symmetries;
    symmetries.push_back(std::make_unique<SpinSymmetry>());
    symmetries.push_back(std::make_unique<InversionSymmetry>());

    experimental::WickOrderedCollector H_normal = experimental::wick_decompose(H, wick_templates);
    experimental::clean_wick_ordered_terms(H_normal, symmetries);
    advanced_clean_up(H_normal);

    // --------------------------------------------------------- //

    const Term cut_generator_half(1, 
        Coefficient::RealInteraction("\\alpha", 
            MomentumList({Momentum('p'), Momentum('q'), Momentum('r')}), 
            IndexWrapper{Index::Sigma, Index::SigmaPrime},
            spin_symmetry
        ),
        SumContainer{MomentumSum{'p', 'q', 'r'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
        std::vector<Operator>({
            c_p_sigma.hermitian_conjugate(),
            c_p_sigma_prime.with_momentum('q').hermitian_conjugate(),
            c_p_sigma_prime.with_momentum(Momentum("q-r")),
            c_p_sigma.with_momentum(Momentum("p+r"))
        })
    );
    Term cut_generator_second_half = cut_generator_half.hermitian_conjugate();
    cut_generator_second_half.flip_sign();

    const std::vector<Term> cut_generator = {cut_generator_half, cut_generator_second_half};

    // --------------------------------------------------------- //

    oss << "We work with the Hamiltonian \n\\begin{align}\nH =" << H 
        << ",\\end{align}\n or in its normal-ordered form\n\\begin{align}H="
        << H_normal << ".\\end{align}\nBased on this, we define the effective single-particle energy\n\\begin{align}"
        << "\\varepsilon(\\mathbf{p}) = \\tilde{\\varepsilon}(\\mathbf{p}) + 2 \\sum_{\\mathbf{q}} \\left[ " 
        << "\\sum_{\\sigma'} U_{\\sigma\\sigma'} (\\mathbf{p}, \\mathbf{q}, 0) \\langle \\hat{n}_{\\mathbf{q},\\sigma'} \\rangle "
        << "- U_{\\sigma \\sigma} (\\mathbf{p}, \\mathbf{q}, \\mathbf{q} - \\mathbf{p}) \\langle \\hat{n}_{\\mathbf{q},\\sigma} \\rangle"
        << " \\right].\n\\end{align}\n"
        << "The generator of the CUT reads\n\\begin{align}\n\\eta =" << cut_generator
        << ".\\end{align}\nHere, we introduced the abbreviation\n\\begin{equation}\n\t"
        << "\\alpha_{\\sigma\\sigma'}(\\mathbf{p}, \\mathbf{q}, \\mathbf{r}) = \\mathrm{sgn} \\left[ "
        << "\\varepsilon(\\mathbf{p}) + \\varepsilon(\\mathbf{q}) - \\varepsilon(\\mathbf{q} - \\mathbf{r}) - \\varepsilon(\\mathbf{p} + \\mathbf{q}) \\right]"
        << " U_{\\sigma \\sigma'} (\\mathbf{p}, \\mathbf{q}, \\mathbf{r})."
        << "\\end{equation}\nThe $\\ell$-dependence of the coefficients is implied. "
        << "Note that, at $\\ell=0$, the interaction is independent of the momenta. " << std::endl;
        
    std::vector<Term> cut_commutator = commutator(H, cut_generator);
    clean_up(cut_commutator);
    
    return cut_commutator;
}

}