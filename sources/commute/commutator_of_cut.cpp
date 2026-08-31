#include "commutator_of_cut.hpp"

#include "../experimental/WickOrderedCollector.hpp"

#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/ExpectationValues>

namespace NickelCUT::commute {
using namespace mrock::symbolic_operators;

const Operator c_p_sigma        = Operator(Momentum('K'), Index::Sigma,      false);
const Operator c_p_sigma_prime  = Operator(Momentum('K'), Index::SigmaPrime, false);

const Operator c_k_up           = Operator(Momentum('k'),  Index::SpinUp,   false);
const Operator c_minus_k_down   = Operator(Momentum("-k"), Index::SpinDown, false);

// It only matters if the indices of the coefficients are identical or differ;
// e.g. U_(sigma, sigma') = U_(simga', sigma), and U_(down, down) = U_(up up)
const std::function<void(Coefficient&)> spin_symmetry = [](Coefficient& coeff) {
    assert(coeff.indices.size() == 2U);
    if (coeff.indices[0] > coeff.indices[1]) {
        std::swap(coeff.indices[0], coeff.indices[1]);
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

TermCollector get_Hamiltonian()
{
    const Term H_kin(1, 
        Coefficient::RealInversionSymmetric("\\tilde{\\varepsilon}", MomentumList(Momentum('K'))),
        SumContainer{MomentumSum{'K'}, IndexSum{ Index::Sigma }},
        std::vector<Operator>({
            c_p_sigma.hermitian_conjugate(),
            c_p_sigma
        })
    );

    const Term H_int(1, 
        Coefficient::RealInteraction("U", MomentumList({Momentum('K'), Momentum('P'), Momentum('Q')}), 
            IndexWrapper{Index::Sigma, Index::SigmaPrime},
            spin_symmetry
        ),
        SumContainer{MomentumSum{'K', 'P', 'Q'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
        std::vector<Operator>({
            c_p_sigma.hermitian_conjugate(),
            c_p_sigma_prime.with_momentum('P').hermitian_conjugate(),
            c_p_sigma_prime.with_momentum(Momentum("P-Q")),
            c_p_sigma.with_momentum(Momentum("K+Q"))
        })
    );

    return {H_kin, H_int};
}

TermCollector get_cut_generator() {
    const Term cut_generator_half(2, 
        Coefficient::RealInteraction("\\alpha", 
            MomentumList({Momentum('K'), Momentum('P'), Momentum('Q')}), 
            IndexWrapper{Index::Sigma, Index::SigmaPrime},
            spin_symmetry
        ),
        SumContainer{MomentumSum{'K', 'P', 'Q'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
        std::vector<Operator>({
            c_p_sigma.hermitian_conjugate(),
            c_p_sigma_prime.with_momentum('P').hermitian_conjugate(),
            c_p_sigma_prime.with_momentum(Momentum("P-Q")),
            c_p_sigma.with_momentum(Momentum("K+Q"))
        })
    );

    //Term cut_generator_second_half = cut_generator_half.hermitian_conjugate();
    //cut_generator_second_half.flip_sign();
    //const TermCollector cut_generator = {cut_generator_half, cut_generator_second_half};

    return {cut_generator_half};
}

TermCollector commutator_of_cut(std::ostringstream& oss) {
    const TermCollector H = get_Hamiltonian();
    const TermCollector cut_generator = get_cut_generator();
    TermCollector cut_commutator = commutator(cut_generator, H);
    cut_commutator.clean_up();

    // The last two terms are identical, but the computer needs special aid to see that
    cut_commutator.back().rename_indices(Index::Sigma, Index::PlaceHolderIndex);
    cut_commutator.back().rename_indices(Index::SigmaPrime, Index::Sigma);
    cut_commutator.back().rename_indices(Index::PlaceHolderIndex, Index::SigmaPrime);
    cut_commutator.back().structure();
    std::sort(cut_commutator.back().sums.spins.begin(), cut_commutator.back().sums.spins.end());

    // Use the spin symmetry of the coefficients
    for(auto& term : cut_commutator) {
        for(auto& coeff : term.coefficients) {
            // The order of indices does not matter
            std::sort(coeff.indices.begin(), coeff.indices.end());
        }
    }
    cut_commutator.combine_duplicates();

    // --------------------------------------------------------- //

    oss << "We work with the Hamiltonian \n\\begin{align}\nH =" << H 
        << ".\\end{align}\nBased on this, we define the effective single-particle energy\n\\begin{align}"
        << "\\varepsilon(\\mathbf{K}) = \\tilde{\\varepsilon}(\\mathbf{K}) + 2 \\sum_{\\mathbf{P}} \\left[ " 
        << "\\sum_{\\sigma'} U_{\\sigma\\sigma'} (\\mathbf{K}, \\mathbf{P}, 0) \\langle \\hat{n}_{\\mathbf{P},\\sigma'} \\rangle "
        << "- U_{\\sigma \\sigma} (\\mathbf{K}, \\mathbf{P}, \\mathbf{P} - \\mathbf{K}) \\langle \\hat{n}_{\\mathbf{P},\\sigma} \\rangle"
        << " \\right].\n\\end{align}\n"
        << "The generator of the CUT reads\n\\begin{align}\n\\eta =" << cut_generator
        << ".\\end{align}\nHere, we introduced the abbreviation\n\\begin{equation}\n\t"
        << "\\alpha_{\\sigma\\sigma'}(\\mathbf{K}, \\mathbf{P}, \\mathbf{Q}) = \\mathrm{sgn} \\left[ "
        << "\\varepsilon(\\mathbf{K}) + \\varepsilon(\\mathbf{P}) - \\varepsilon(\\mathbf{P} - \\mathbf{Q}) - \\varepsilon(\\mathbf{K} + \\mathbf{P}) \\right]"
        << " U_{\\sigma \\sigma'} (\\mathbf{K}, \\mathbf{P}, \\mathbf{Q})."
        << "\\end{equation}\nThe $\\ell$-dependence of the coefficients is implied. "
        << "Note that, at $\\ell=0$, the interaction is independent of the momenta. " << std::endl;
    
    oss << "\nThe commutator reads\n\\begin{align}\n[\\eta, H] =" 
        //<< "\\text{1 page of terms}" 
        << cut_commutator
        << "\\end{align}" << std::endl;

    return cut_commutator;
}

TermCollector commutator_pair_annihilation(std::ostringstream& oss) {
    const Term pair_annihilation(1, 
        std::vector<Operator>({
            c_minus_k_down,
            c_k_up
        })
    );

    const TermCollector cut_generator = get_cut_generator();
    TermCollector cut_commutator = commutator(cut_generator, pair_annihilation);
    cut_commutator.clean_up();

    oss << "Next, we watch how the pair annihilation (and by extend the pair creation operators) transform due to the flow."
        << "To do so, we evaluate\n\\begin{align*}"
        << "[\\eta (\\ell), O (\\ell)] =" << cut_commutator
        << "\\end{align*}\n, where $O(\\ell)$ represents the pair annihilation operator at $\\ell$.";

    return cut_commutator;
}

}