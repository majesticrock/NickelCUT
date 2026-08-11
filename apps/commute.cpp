#include "../sources/experimental/WickOrderedCollector.hpp"
#include "../sources/commute/commutator_of_cut.hpp"
#include "../sources/commute/extract_flow_coefficients.hpp"
#include "../sources/commute/restructure_terms.hpp"
#include "../sources/commute/export_as_flow_equation.hpp"

#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/Wick.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>
#include <mrock/symbolic_operators/WickSymmetry.hpp>
#include <mrock/symbolic_operators/SerializationHeaders.hpp>

#include <mrock/utility/OutputConvenience.hpp>

#include <cstddef>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <memory>
#include <vector>
#include <list>
#include <iostream>

using namespace mrock::symbolic_operators;
using namespace NickelCUT::commute;

int main(){
    std::ostringstream oss;

    const std::vector<Term> commutator = commutator_of_cut(oss);
    oss << "\nThe commutator reads\n\\begin{align}\n[H, \\eta] =" 
        << "\\text{2 pages of terms}" 
        // << commutator
        << "\\end{align}" << std::endl;

    const std::vector<WickOperatorTemplate> wick_templates({ 
        WickOperatorTemplate({Num_Comparison}, Momentum(), OperatorType::Number) 
    });

    std::vector<std::unique_ptr<WickSymmetry>> symmetries;
    symmetries.push_back(std::make_unique<SpinSymmetry>());
    symmetries.push_back(std::make_unique<InversionSymmetry>());


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
        << "\\partial_\\ell \\varepsilon (\\mathbf{p}) ="
        << flow_coefficients[0] 
        << "\\end{align*}" << std::endl;
    
    oss << "\\begin{align*}\n\t"
        << "\\partial_\\ell U_{\\uparrow \\downarrow} (\\mathbf{q}, \\mathbf{p}, \\mathbf{r}) ="
        << flow_coefficients[1] 
        << "\\end{align*}" << std::endl;

    oss << "\\begin{align*}\n\t"
        << "\\partial_\\ell U_{\\parallel} (\\mathbf{q}, \\mathbf{p}, \\mathbf{r}) ="
        << flow_coefficients[2] 
        << "\\end{align*}" << std::endl;

    mrock::utility::save_string_raw(oss.str(), "commute_output.tex");

    export_as_flow_equation(flow_coefficients);

    return 0;
}