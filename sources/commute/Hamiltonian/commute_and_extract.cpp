#include "commute_and_extract.hpp"

#include "extract_flow_coefficients.hpp"
#include "restructure_terms.hpp"
#include "export_as_flow_equation.hpp"

#include "../verify/Verifier.hpp"
#include "../commutator_of_cut.hpp"
#include "../../experimental/WickOrderedCollector.hpp"

#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/ExpectationValues>

#include <vector>
#include <list>

//#define RUN_VERFICATION

namespace NickelCUT::commute::Hamiltonian
{
using namespace mrock::symbolic_operators;
using namespace NickelCUT::commute;

experimental::WickOrderedCollector commute_and_normal_order(std::ostringstream& oss)
{
    const TermCollector commutator = commutator_of_cut(oss);

    const std::vector<WickOperatorTemplate> wick_templates = get_wick_templates();
    const std::vector<std::unique_ptr<WickSymmetry>> symmetries = get_symmetries();

    experimental::WickOrderedCollector normal_ordered_result = experimental::wick_decompose(commutator, wick_templates);
    experimental::clean_wick_ordered_terms(normal_ordered_result, symmetries);

#ifdef RUN_VERFICATION
    const verify::Verifier verifier;
    verify::Verifier::SparseMatrix matrix_rep = verifier.symbolic_to_matrix(commutator);

    if (!verifier.matrices_equal(matrix_rep, matrix_rep.transpose())) {
        throw std::runtime_error("Matrix representation is not Hermitian!");
    }
    else {
        std::cout << "Initial commutator is Hermitian." << std::endl;
    }

    if (!verifier(normal_ordered_result, "CUT commutator after basic cleaning", matrix_rep)) {
        throw std::runtime_error("Matrix representations do not match!");
    }
    else {
        std::cout << "Basic cleaning is okay." << std::endl;
    }
#endif

    // We don't care about sextic terms (delete size > 4) and we also don't care about the
    // renormalization of the groundstate energy (delte empty expressions)
    std::erase_if(normal_ordered_result.terms, [](const experimental::WickOrderedTerm& term){
        return term.wick_expression.size() > 4U || term.wick_expression.empty();
    });

#ifdef RUN_VERFICATION
    // The matrix must change as we delete all sextic and higher terms.
    matrix_rep = verifier.symbolic_to_matrix(normal_ordered_result);

    advanced_clean_up(normal_ordered_result);

    if (!verifier(normal_ordered_result, "CUT commutator after advanced cleaning", matrix_rep)) {
        throw std::runtime_error("Matrix representations do not match!");
    }
    else {
        std::cout << "Advanced clean up is okay." << std::endl;
    }
#endif

    oss << "After normal ordering with respect to the Fermi sea, we omit any contribution with more than 4 operators. "
        << "The result reads\n\\begin{align*}\n\t"
        << "\\text{3 pages of terms}"
        //<< normal_ordered_result 
        << "\\end{align*}" << std::endl;
    
    return normal_ordered_result;
}

void commute_and_extract(std::ostringstream& oss) {
    experimental::WickOrderedCollector normal_ordered_result = commute_and_normal_order(oss); 
    oss << "The unique types of Wick-ordered expressions are\n\\begin{align*}\n";
    std::list<experimental::WickOrderedExpression> unique_wicks;
    for(const auto& term : normal_ordered_result) {
        if (!exists_in(unique_wicks, term.wick_expression)) {
            unique_wicks.push_back(term.wick_expression);
            oss << "\t&" << term.wick_expression << "\\\\\n";
        }
    }
    oss << "\\end{align*}" << std::endl;

    std::array<WickTermCollector, 3> flow_coefficients = extract_flow_coefficients(normal_ordered_result);

//#ifdef RUN_VERIFICATION
    const verify::Verifier verifier;
    verify::Verifier::SparseMatrix compare_matrix = verifier.symbolic_to_matrix(normal_ordered_result);
    std::cout << "Cached compare_matrix." << std::endl;

    // In the next step we drop all the normal ordered terms because only one structure of terms remains
    // For comparing them to our matrix result, however, we need to readd those.
    experimental::WickOrderedCollector readded;

    auto readd_wick_expressions = [&readded](const WickTerm& term, int n_coeff) {
        readded.push_back(experimental::WickOrderedTerm(term));
  
        readded.back().sums.momenta.push_back('K');
        if (n_coeff == 0 || n_coeff == 2) {
            readded.back().sums.spins.push_back(Index::Sigma);
        }
        if (n_coeff == 1 || n_coeff == 2) {
            readded.back().sums.momenta.push_back('P');
            readded.back().sums.momenta.push_back('Q');
        }

        if (n_coeff == 0) {
            readded.back().wick_expression.operators = {
                Operator(Momentum('K'), Index::Sigma, true),
                Operator(Momentum('K'), Index::Sigma, false),
            };
        }
        else if (n_coeff == 1) {
            readded.back().wick_expression.operators = {
                Operator(Momentum('K'), Index::SpinUp, true),
                Operator(Momentum('P'), Index::SpinDown, true),
                Operator(Momentum("P-Q"), Index::SpinDown, false),
                Operator(Momentum("K+Q"), Index::SpinUp, false)
            };

            readded.push_back(readded.back());

            readded.back().wick_expression[0].indices[0] = Index::SpinDown;
            readded.back().wick_expression[1].indices[0] = Index::SpinUp;
            readded.back().wick_expression[2].indices[0] = Index::SpinUp;
            readded.back().wick_expression[3].indices[0] = Index::SpinDown;
        }
        else if (n_coeff == 2) {
            readded.back().wick_expression.operators = {
                Operator(Momentum('K'), Index::Sigma, true),
                Operator(Momentum('P'), Index::Sigma, true),
                Operator(Momentum("P-Q"), Index::Sigma, false),
                Operator(Momentum("K+Q"), Index::Sigma, false)
            };
        }
    };

    for(int i=0; i<3; ++i) {
        for (const auto& term : flow_coefficients[i]) {
            readd_wick_expressions(term, i);
        }
    }
    if(!verifier(readded, "Extracted flow coefficients", compare_matrix)) {
        throw std::runtime_error("Matrix representations do not match!");
    }
    else {
        std::cout << "Extracted flow coefficients are okay." << std::endl;
    }
//#endif

    for (auto& fc : flow_coefficients) {
        improve_flow_coefficient_structure(fc);
    }

//#ifdef
    readded.clear();
    for(int i=0; i<3; ++i) {
        for (const auto& term : flow_coefficients[i]) {
            readd_wick_expressions(term, i);
        }
    }
    if(!verifier(readded, "Extracted flow coefficients after restructuring", compare_matrix)) {
        throw std::runtime_error("Matrix representations do not match!");
    }
    else {
        std::cout << "Extracted flow coefficients after restructuring are okay." << std::endl;
    }
//#endif

    oss << "This leaves the flow of the coefficients as follows:\n\\begin{align*}\n\t"
        << "\\partial_\\ell \\varepsilon (\\mathbf{K}) ="
        << flow_coefficients[0] 
        << "\\end{align*}" << std::endl;
    
    oss << "\\begin{align*}\n\t"
        << "\\partial_\\ell U_{\\uparrow \\downarrow} (\\mathbf{K}, \\mathbf{P}, \\mathbf{Q}) ="
        << flow_coefficients[1] 
        << "\\end{align*}" << std::endl;

    oss << "\\begin{align*}\n\t"
        << "\\partial_\\ell U_{\\parallel} (\\mathbf{K}, \\mathbf{P}, \\mathbf{Q}) ="
        << flow_coefficients[2] 
        << "\\end{align*}" << std::endl;

    //////////////////////////////////////////////////////////////
    export_as_flow_equation(flow_coefficients);
}
} // namespace NickelCUT::commute::Hamiltonian


