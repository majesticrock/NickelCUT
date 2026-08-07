#include "../sources/experimental/WickOrderedCollector.hpp"
#include "../sources/commute_helper.hpp"

#include <mrock/symbolic_operators/Commutation>
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
    std::cout << "\nThe commutator reads\n\\begin{align}\n[H, \\eta] =" << commutator
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

    for (auto& term : normal_ordered_result) {
        if (term.wick_expression.empty()) {
            restructure_expectation_values(term, Index::Sigma, {}, 'p', {});
        }
        else if (term.is_bilinear()) {
            restructure_bilinear_term(term);
        }
        else if (term.is_quartic()) {
            restructure_quartic_term(term);
        }

        for (auto& coeff : term.coefficients) {
            coeff.use_inversion_symmetry();
        }
    }

    normal_ordered_result.combine_duplicates();

    std::cout << "After normal ordering with respect to the Fermi sea, we omit any contribution with more than 4 operators. "
        << "The result reads\n\\begin{align*}\n\t"
        << normal_ordered_result << "\\end{align*}" << std::endl;

    std::cout << "The unique types of Wick-ordered expressions are\n\\begin{align*}\n";

    std::list<experimental::WickOrderedExpression> unique_wicks;
    for(const auto& term : normal_ordered_result) {
        if (!exists_in(unique_wicks, term.wick_expression)) {
            unique_wicks.push_back(term.wick_expression);
            std::cout << "\t&" << term.wick_expression << "\\\\\n";
        }
    }
    std::cout << "\\end{align*}" << std::endl;

    std::array<WickTerm, 3> flow_coefficients;

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

    return 0;
}