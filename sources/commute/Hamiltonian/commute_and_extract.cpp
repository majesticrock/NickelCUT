#include "commute_and_extract.hpp"

#include "extract_flow_coefficients.hpp"
#include "restructure_terms.hpp"
#include "export_as_flow_equation.hpp"

#include "../commutator_of_cut.hpp"
#include "../../experimental/WickOrderedCollector.hpp"

#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/ExpectationValues>

#include <vector>
#include <list>

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

    std::erase_if(normal_ordered_result.terms, [](const experimental::WickOrderedTerm& term){
        return term.wick_expression.size() > 4U;
    });

    advanced_clean_up(normal_ordered_result);

    oss << "After normal ordering with respect to the Fermi sea, we omit any contribution with more than 4 operators. "
        << "The result reads\n\\begin{align*}\n\t"
        << "\\text{10 pages of terms}"
        // << normal_ordered_result 
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
    for (auto& fc : flow_coefficients)
        improve_flow_coefficient_structure(fc);

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


