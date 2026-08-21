#include "../sources/commute/commutator_of_cut.hpp"
#include "../sources/commute/Hamiltonian/commute_and_extract.hpp"
#include "../sources/commute/Hamiltonian/restructure_terms.hpp"
#include "../sources/experimental/WickOrderedCollector.hpp"

#include <mrock/utility/OutputConvenience.hpp>

#include <sstream>

using namespace NickelCUT::commute;
using namespace mrock::symbolic_operators;
using namespace mrock::symbolic_operators::experimental;

int main(int /*argc*/, char** /*argv*/) {
    std::ostringstream oss;

    const std::vector<Term> H = get_Hamiltonian();
    const std::vector<Term> cut_generator = get_cut_generator();

    const std::vector<WickOperatorTemplate> wick_templates = get_wick_templates();
    const std::vector<std::unique_ptr<WickSymmetry>> symmetries = get_symmetries();

    WickOrderedCollector H_norm = wick_decompose(H, wick_templates);
    clean_wick_ordered_terms(H_norm, symmetries);
    Hamiltonian::advanced_clean_up(H_norm);

    WickOrderedCollector H_diff = H_norm - H_norm.hermitian_conjugate();
    Hamiltonian::advanced_clean_up(H_diff);

    oss << "We see that\\begin{align*}\n\tH - H^\\dagger =" 
        << H_diff 
        << ".\\end{align*}\n";

    WickOrderedCollector eta_norm = wick_decompose(cut_generator, wick_templates);
    clean_wick_ordered_terms(eta_norm, symmetries);
    Hamiltonian::advanced_clean_up(eta_norm);

    WickOrderedCollector eta_diff = eta_norm - eta_norm.hermitian_conjugate();
    Hamiltonian::advanced_clean_up(eta_diff);
    oss << "Similar, we find\\begin{align*}\n\t\\eta + \\eta^\\dagger =" 
        << eta_diff 
        << ".\\end{align*}\n";

    //
    //const WickOrderedCollector commutator_eta_H = Hamiltonian::commute_and_normal_order(oss);
    //
    ////////////////////////////////////////////////////////////////
    //WickOrderedCollector hermitian_difference = commutator_eta_H - commutator_eta_H.hermitian_conjugate();
    //Hamiltonian::advanced_clean_up(hermitian_difference);

    
    
    //oss << "Moreover, we find\n\\begin{align*}\n\t[H, \\eta] - [H, \\eta]^\\dagger =" 
    //    << hermitian_difference 
    //    << ".\\end{align*}\n";

    mrock::utility::save_string_raw(oss.str(), "test_output.tex");

    return 0;
}