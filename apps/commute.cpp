#include "../sources/experimental/WickOrderedCollector.hpp"
#include "../sources/commute/commutator_of_cut.hpp"
#include "../sources/commute/extract_flow_coefficients.hpp"
#include "../sources/commute/restructure_terms.hpp"

#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/Wick.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>
#include <mrock/symbolic_operators/WickSymmetry.hpp>
#include <mrock/symbolic_operators/SerializationHeaders.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <vector>
#include <list>
#include <iostream>

using namespace mrock::symbolic_operators;
using namespace NickelCUT::commute;

int main(){
    const std::vector<Term> commutator = commutator_of_cut();
    std::cout << "\nThe commutator reads\n\\begin{align}\n[H, \\eta] =" 
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

    std::cout << "After normal ordering with respect to the Fermi sea, we omit any contribution with more than 4 operators. "
        << "The result reads\n\\begin{align*}\n\t"
        << "\\text{10 pages of terms}"
        // << normal_ordered_result 
        << "\\end{align*}" << std::endl;

    std::cout << "The unique types of Wick-ordered expressions are\n\\begin{align*}\n";

    std::list<experimental::WickOrderedExpression> unique_wicks;
    for(const auto& term : normal_ordered_result) {
        if (!exists_in(unique_wicks, term.wick_expression)) {
            unique_wicks.push_back(term.wick_expression);
            std::cout << "\t&" << term.wick_expression << "\\\\\n";
        }
    }
    std::cout << "\\end{align*}" << std::endl;

    std::array<WickTermCollector, 3> flow_coefficients = extract_flow_coefficients(normal_ordered_result);
    for (auto& fc : flow_coefficients)
        improve_flow_coefficient_structure(fc);

    // serialization
    std::filesystem::create_directories("apps/commute_output/");
    // create an output file stream and a text archive to serialize the vector
    {
        std::ofstream ofs("apps/commute_output/identity.bin", std::ios::binary);
        boost::archive::binary_oarchive oa(ofs);
        oa << flow_coefficients[0];
        ofs.close();
    }
    {
        std::ofstream ofs("apps/commute_output/bilinear.bin", std::ios::binary);
        boost::archive::binary_oarchive oa(ofs);
        oa << flow_coefficients[1];
        ofs.close();
    }
    {
        std::ofstream ofs("apps/commute_output/quartic.bin", std::ios::binary);
        boost::archive::binary_oarchive oa(ofs);
        oa << flow_coefficients[2];
        ofs.close();
    }

    std::cout << "This leaves the flow of the coefficients as follows:\n\\begin{align*}\n\t"
        << "\\partial_\\ell C_0 ="
        << "\\text{3 pages of terms}"
        //<< flow_coefficients[0] 
        << "\\end{align*}" << std::endl;
    
    std::cout << "\\begin{align*}\n\t"
        << "\\partial_\\ell \\varepsilon (\\mathbf{p}) ="
        << flow_coefficients[1] 
        << "\\end{align*}" << std::endl;

    std::cout << "\\begin{align*}\n\t"
        << "\\partial_\\ell U_{\\sigma, \\sigma'} (\\mathbf{q}, \\mathbf{p}, \\mathbf{r}) ="
        << flow_coefficients[2] 
        << "\\end{align*}" << std::endl;

    return 0;
}