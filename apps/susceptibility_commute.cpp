#include "../sources/commute/commutator_of_cut.hpp"

#include <mrock/symbolic_operators/Commutation>
#include <mrock/symbolic_operators/ExpectationValues>

#include <iostream>
#include <vector>

using namespace mrock::symbolic_operators;

const Momentum base_k = Momentum('k');
const Momentum base_x = Momentum('x');

const Operator c_k = Operator{base_k, Index::SpinUp, false};
const Operator c_minus_k = Operator{-base_k, Index::SpinDown, false};
const Operator c_k_dagger = Operator{base_k, Index::SpinUp, true};
const Operator c_minus_k_dagger = Operator{-base_k, Index::SpinDown, true};

std::vector<TermCollector> get_basis_operators() 
{
    std::vector<TermCollector> ret = {// 0/1: f, f^+
            TermCollector({Term(1, std::vector<Operator>({c_minus_k, c_k}))}),
            TermCollector({Term(1, std::vector<Operator>({c_k_dagger, c_minus_k_dagger}))}),
            // 2/3: n_up/down
            TermCollector({Term(1, std::vector<Operator>({c_k_dagger, c_k}))}),
            TermCollector({Term(1, std::vector<Operator>({c_minus_k_dagger, c_minus_k}))}),
    };

    //for (auto& _v : ret) {
    //    for (auto& v : _v) {
    //        if (v.operators.front().is_daggered) {
    //            v.operators.front().momentum += Momentum('x');
    //        } else {
    //            v.operators.front().momentum += Momentum('x', -1);
    //        }
    //    }
    //}
    return ret;
}

void carry_out_spin_summations(WickTermCollector& terms) {
    const std::size_t size = terms.size();
    for (std::size_t i=0U; i < size; ++i) {
        if (terms[i].sums.spins.empty()) continue;
        if (terms[i].sums.spins.size() > 1U) throw std::invalid_argument("Expected no more than 1 spin summation!");

        terms[i].sums.spins.clear();
        terms.push_back(terms[i]);
        
        terms[i].replace_each_index(Index::Sigma, Index::SpinDown);
        terms.back().replace_each_index(Index::Sigma, Index::SpinUp);
    }

    for (auto& term : terms) {
        IndexWrapper& coeff_indices = term.coefficients.front().indices;
        if (coeff_indices.size() != 2U) continue;

        if (coeff_indices[0] == coeff_indices[1]) {
            coeff_indices.clear();
            coeff_indices.push_back(Index::Parallel);
        }
        else {
            coeff_indices.clear();
            coeff_indices.push_back(Index::AntiParallel);
        }

        term.coefficients.front().custom_symmetry = std::nullopt;
    }
}

int main(int, char**)
{
    const auto H = NickelCUT::commute::get_Hamiltonian();
    const auto wick_templates = {WickOperatorTemplate{{SC_Comparison}, Momentum(), OperatorType::SC},
            WickOperatorTemplate{{SC_Comparison}, Momentum(std::vector<MomentumSymbol>(), true), OperatorType::Eta},
            WickOperatorTemplate{{Num_Comparison}, Momentum(), OperatorType::Number},
            WickOperatorTemplate{{Num_Comparison}, Momentum(std::vector<MomentumSymbol>(), true), OperatorType::CDW}};
    const auto symmetries = NickelCUT::commute::get_symmetries();

    const std::vector<TermCollector> basis = get_basis_operators();
    std::vector<TermCollector> basis_daggered(basis);
    for (auto& t : basis_daggered) {
        t.hermitian_conjugate();
        t.rename_momenta('k', 'l');
    }

    std::cout << "\\begin{align*}\n\t H =" << H << "\\end{align*}\n" << std::endl;

    for (std::size_t i = 0U; i < basis.size(); ++i) {
        TermCollector commute_with_H = commutator(H, basis[i]);
        commute_with_H.clean_up();

        for (std::size_t j = 0U; j < basis.size(); ++j) {
            std::cout << "\\subsection{" << i << "." << j << "}" << std::endl;
            TermCollector terms = commutator(basis_daggered[j], commute_with_H);
            terms.clean_up();

            WickTermCollector wicks;
            wicks_theorem(terms, wick_templates, wicks);
            wicks.clean_up(symmetries);
            //carry_out_spin_summations(wicks);
            //wicks.clean_up(symmetries);

            std::cout << "\\begin{align*}\n\t\\langle [ " << basis_daggered[j].to_string_without_prefactor()
                        << ", [H, " << basis[i].to_string_without_prefactor() << " ]] \\rangle =" << wicks
                        << "\\end{align*}" << std::endl;

            terms.clear();
            wicks.clear();
            terms = commutator(basis_daggered[j], basis[i]);
            terms.clean_up();
            wicks_theorem(terms, wick_templates, wicks);
            //wicks.clean_up(symmetries);

            std::cout << "\\begin{align*}\n\t[ " << basis_daggered[j].to_string_without_prefactor() << ", "
                        << basis[i].to_string_without_prefactor() << " ] =" << terms << "\\end{align*}" << std::endl;
            std::cout << "\\begin{align*}\n\t[ " << basis_daggered[j].to_string_without_prefactor() << ", "
                        << basis[i].to_string_without_prefactor() << " ] =" << wicks << "\\end{align*}" << std::endl;
        }
    }

    return 0;
}