#include "commutator_of_cut.hpp"

#include "../experimental/WickOrderedCollector.hpp"

#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/ExpectationValues>

namespace NickelCUT::commute {
using namespace mrock::symbolic_operators;

const Operator c_p_sigma        = Operator(Momentum('p'), Index::Sigma,      false);
const Operator c_p_sigma_prime  = Operator(Momentum('p'), Index::SigmaPrime, false);

const Operator c_k_up           = Operator(Momentum('k'),  Index::SpinUp,   false);
const Operator c_minus_k_down   = Operator(Momentum("-k"), Index::SpinDown, false);

// It only matters if the indizes of the coefficients are identical or differ;
// e.g. U_(sigma, sigma') = U_(simga', sigma), and U_(down, down) = U_(up up)
const std::function<void(Coefficient&)> spin_symmetry = [](Coefficient& coeff) {
    assert(coeff.indizes.size() == 2U);
    if (coeff.indizes[0] > coeff.indizes[1]) {
        std::swap(coeff.indizes[0], coeff.indizes[1]);
    }
};

std::vector<WickOperatorTemplate> get_wick_templates()
{ 
    return { WickOperatorTemplate({Num_Comparison}, Momentum(), OperatorType::Number) };
}

std::vector<std::unique_ptr<WickSymmetry>> get_symmetries() {
    std::vector<std::unique_ptr<WickSymmetry>> symmetries;
    symmetries.push_back(std::make_unique<SpinSymmetry>());
    symmetries.push_back(std::make_unique<InversionSymmetry>());
    return symmetries;
}

std::vector<Term> get_Hamiltonian()
{
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

    return {H_kin, H_int};
}

std::vector<Term> get_cut_generator() {
    const Term cut_generator_half(2, 
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

    //Term cut_generator_second_half = cut_generator_half.hermitian_conjugate();
    //cut_generator_second_half.flip_sign();
    //const std::vector<Term> cut_generator = {cut_generator_half, cut_generator_second_half};

    return {cut_generator_half};
}

std::vector<Term> commutator_of_cut(std::ostringstream& oss) {
    const std::vector<Term> H = get_Hamiltonian();
    const std::vector<Term> cut_generator = get_cut_generator();
    std::vector<Term> cut_commutator = commutator(cut_generator, H);
    clean_up(cut_commutator);

    // --------------------------------------------------------- //

    oss << "We work with the Hamiltonian \n\\begin{align}\nH =" << H 
        << ".\\end{align}\nBased on this, we define the effective single-particle energy\n\\begin{align}"
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
    
    oss << "\nThe commutator reads\n\\begin{align}\n[\\eta, H] =" 
        //<< "\\text{2 pages of terms}" 
        << cut_commutator
        << "\\end{align}" << std::endl;

    return cut_commutator;
}

std::vector<Term> commutator_pair_annihilation(std::ostringstream& oss) {
    const Term pair_annihilation(1, 
        std::vector<Operator>({
            c_minus_k_down,
            c_k_up
        })
    );

    const std::vector<Term> cut_generator = get_cut_generator();
    std::vector<Term> cut_commutator = commutator(cut_generator, pair_annihilation);
    clean_up(cut_commutator);

    oss << "Next, we watch how the pair annihilation (and by extend the pair creation operators) transform due to the flow."
        << "To do so, we evaluate\n\\begin{align*}"
        << "[\\eta (\\ell), O (\\ell)] =" << cut_commutator
        << "\\end{align*}\n, where $O(\\ell)$ represents the pair annihilation operator at $\\ell$.";

    return cut_commutator;
}

}