#include "export_as_flow_equation.hpp"
#include "../../experimental/WickOrderedCollector.hpp"

#include <mrock/utility/OutputConvenience.hpp>

#include <algorithm>
#include <cmath>
#include <cassert>
#include <string>
#include <iostream>

namespace NickelCUT::commute::Hamiltonian {
using namespace mrock::symbolic_operators;

const std::string outer_K_loop = "#pragma omp parallel for\nfor (int K_pos=0; K_pos < N; ++K_pos) {\nmomentum_iterator<L> K(K_pos);\n";

std::string momentum_for_loop(std::string const& it_name) {
    return "for (momentum_iterator<L> " + it_name + " = momentum_iterator<L>::begin(); " 
        + it_name + " != momentum_iterator<L>::end(); ++" + it_name + ") {\n";
}

std::string sign_to_string(int factor) {
    return (factor < 0 ? "-" : "+");
}

std::string momentum_to_code(Momentum momentum, bool make_positive=false) {
    if (momentum.is_zero()) {
        return "Gamma<L>";
    }

    if (make_positive && momentum.first_momentum_is_negative()) {
        momentum.flip_momentum();
    }
    std::string code = "";
    for (const auto& symbol : momentum) {
        assert(symbol.factor == 1 || symbol.factor == -1);
        if (code == "") {
            if (symbol.factor == -1) code += "-";
        }
        else {
            code += sign_to_string(symbol.factor);
        }
        code += symbol.name._n;
    }
    return code;
}

std::string access_coefficient(const Coefficient& coeff) {
    std::string code = "";
    if (coeff.name == "\\alpha") {
        //code += "sign(";
        //
        //code += "current.epsilon_tilde[";
        //code += momentum_to_code(coeff.momenta[0], true);
        //code += "] + ";
        //code += "current.epsilon_tilde[";
        //code += momentum_to_code(coeff.momenta[1], true);
        //code += "] - ";
//
        //code += "current.epsilon_tilde[";
        //code += momentum_to_code(coeff.momenta[0] + coeff.momenta[2], true);
        //code += "] - ";
        //code += "current.epsilon_tilde[";
        //code += momentum_to_code(coeff.momenta[1] - coeff.momenta[2], true);
        //code += "]";
//
        //code += ") \n\t* ";

        code += "alpha_sign_cache(";
        code += momentum_to_code(coeff.momenta[0]) + ", " + momentum_to_code(coeff.momenta[1]) + ", " + momentum_to_code(coeff.momenta[2]);
        code += ") \n\t* ";
    }
    if (coeff.name == "U" || coeff.name == "\\alpha") {
        if (coeff.indices[0] == Index::Parallel) {
            code += "current.interactions_same_spin";
        }
        else {
            code += "current.interactions_differing_spin";
        }
        code += "(";
        for (auto it = coeff.momenta.begin(); it != coeff.momenta.end(); ++it) {
            code += momentum_to_code(*it);
            if (it != coeff.momenta.end() - 1) code += ", ";
        }
        code += ")";
    }
    if (coeff.name == "\\tilde{\\varepsilon}") {
        code += "current.epsilon_tilde[";
        code += momentum_to_code(coeff.momenta[0]);
        code += "]";
    }

    return code;
}

IntFractional get_smallest_factor(const experimental::WickOrderedCollector& terms) {
    IntFractional smallest_factor = std::min_element(terms.begin(), terms.end(), [](const auto& l, const auto& r) {
        return std::abs((double)l.multiplicity) < std::abs((double)r.multiplicity);
    })->multiplicity;
    if (smallest_factor.numerator < 0) smallest_factor.numerator *= -1;
    return smallest_factor;
}

std::string generate_bilinear(const experimental::WickOrderedCollector& bilinears) 
{
    const std::string accessor = "dHdl.dispersion[K]";

    std::string code = outer_K_loop;
    code += momentum_for_loop("P");
    code += "double nQ_value{};\ndouble one_value{};\n";
    code += momentum_for_loop("Q");

    const IntFractional smallest_factor = get_smallest_factor(bilinears);

    for (const auto& term : bilinears) {
        assert(term.operators[0].momentum == Momentum('P'));
        assert(term.coefficients[0].name == "U");
        assert(term.coefficients[1].name == "\\alpha");

        if (term.operators.size() == 1U) {
            code += "one_value ";
        }
        else {
            assert(term.operators[1].momentum == Momentum('Q'));
            code += "nQ_value ";
        }

        const auto mult = term.multiplicity / smallest_factor;
        code += mult > 0 ? "+= " : "-= ";
        if (mult != 1 && mult != -1) {
            code += std::to_string(std::abs(static_cast<double>(mult)));
            code += " * ";
        }
        
        code += access_coefficient(term.coefficients[0]);
        code += "\n\t* ";
        code += access_coefficient(term.coefficients[1]);
        code += ";\n";
    }

    code += "nQ_value *= occupation_numbers[Q];\n";
    code += "} // Q-loop\n";
    code += accessor + " += " + std::to_string(static_cast<double>(smallest_factor)) 
        + " * (nQ_value + one_value) * occupation_numbers[P];\n";
    code += "} // P-loop\n";
    code += "} // K-loop\n";
    return code;
}

std::string generate_quartic(const experimental::WickOrderedCollector& quartics, bool parallel) {
    const std::string accessor = std::string("dHdl.") 
        + (parallel ? std::string("interactions_same_spin") : std::string("interactions_differing_spin"))
        + std::string("(K, P, Q)");

    std::string code = outer_K_loop;
    code += momentum_for_loop("P");
    if (parallel) {
        code += "if (K==P) continue; // Pauli principle\n";
    }
    code += momentum_for_loop("Q");
    if (parallel) {
        code += "if (K+Q==P-Q) continue; // Pauli principle\n";
    }

    for (auto it = quartics.begin(); it != quartics.end() && it->sums.momenta.empty(); ++it) {
        code += accessor + (it->multiplicity > 0 ? "+= " : "-= ");
        if (it->multiplicity != 1 && it->multiplicity != -1) {
            code += std::to_string(std::abs(static_cast<double>(it->multiplicity)));
            code += " * ";
        }

        code += access_coefficient(it->coefficients[0]);
        code += "\n\t* ";
        code += access_coefficient(it->coefficients[1]);

        if (!(it->operators.empty())) {
            assert(it->operators.size() == 1U);
            code += "\n\t* occupation_numbers[";
            code += momentum_to_code(it->operators[0].momentum);
            code += "]";
        }
        code += ";\n";
    }

    code += momentum_for_loop("R");
    code += "double nR_value{};\ndouble one_value{};\n";

    const IntFractional smallest_factor = get_smallest_factor(quartics);

    for (auto& term : quartics) {
        if (term.sums.momenta.empty()) continue;

        code += term.operators.empty() ? "one_value " : "nR_value ";

        const auto mult = term.multiplicity / smallest_factor;
        code += mult > 0 ? "+= " : "-= ";
        if (mult != 1 && mult != -1) {
            code += std::to_string(std::abs(static_cast<double>(mult)));
            code += " * ";
        }

        code += access_coefficient(term.coefficients[0]);
        code += "\n\t* ";
        code += access_coefficient(term.coefficients[1]);

        if (term.operators.size() == 2U) {
            code += "\n\t* occupation_numbers[";
            if (term.operators[0].momentum == Momentum('R')) {
                code += momentum_to_code(term.operators[1].momentum);
            }
            else if (term.operators[1].momentum == Momentum('R')) {
                code += momentum_to_code(term.operators[0].momentum);
            }
            else {
                code += "QUARTIC ONE SHOULD BE R!";
                return code;
            }
            code += "]";
        }
        code += ";\n";
    }

    code += accessor + " += " + std::to_string(static_cast<double>(smallest_factor)) 
        + " * (one_value + occupation_numbers[R] * nR_value);\n";

    code += "} // R-loop\n";

    code += "} // Q-loop\n";
    code += "} // P-loop\n";
    code += "} // K-loop\n";

    return code;
}

void export_as_flow_equation(const std::array<experimental::WickOrderedCollector, 3> flow_coeffs) 
{
    const std::string file_header = 
        "#include \"FlowEquation.hpp\"\n\n"
        "#include \"momentum_iterator.hpp\"\n"
        "#include \"FlowContainer.hpp\"\n"
        "#include \"occupation_numbers.hpp\"\n"
        "#include \"../helper_functions.hpp\"\n"
        "\n"
        "#include <omp.h>\n"
        "\n"
        "namespace NickelCUT::flow {\n\n"
        "void FlowEquation::operator()(const FlowContainer& current, FlowContainer& dHdl, const double /*l*/) {\n"
        "dHdl.reset();\n"
        "static InteractionDataFrame alpha_sign_cache;\n"
        "#pragma omp parallel for\n"
        "for (int K_pos=0; K_pos < N; ++K_pos) {\n"
        "momentum_iterator<L> K(K_pos);\n"
        "for (momentum_iterator<L> P = momentum_iterator<L>::begin(); P != momentum_iterator<L>::end(); ++P) {\n"
        "for (momentum_iterator<L> Q = momentum_iterator<L>::begin(); Q != momentum_iterator<L>::end(); ++Q) {\n"
        "alpha_sign_cache(K,P,Q) = sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[P-Q] - current.epsilon_tilde[K+Q]);\n"
        "} // Q-loop\n"
        "} // P-loop\n"
        "} // K-loop\n\n";

    const std::string file_footer = 
        "dHdl.interactions_same_spin.antisymmetrize();\n"
        "dHdl.interactions_same_spin.symmetrize();\n"
        "dHdl.interactions_differing_spin.symmetrize();\n"
        "}\n} // namespace NickelCUT::flow";

    const std::string bilinear_code = generate_bilinear(flow_coeffs[0]) 
        + "\n//--------------------------------------------------------------//\n";
    const std::string u_anti_code = generate_quartic(flow_coeffs[1], false)
        + "\n//--------------------------------------------------------------//\n";
    const std::string u_para_code = generate_quartic(flow_coeffs[2], true);

    mrock::utility::save_string_raw(file_header + bilinear_code + u_anti_code + u_para_code + file_footer, "sources/flow/FlowEquation.cpp");
}

}