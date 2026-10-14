#include "export_as_flow_equation.hpp"
#include <mrock/symbolic_operators/experimental/WickOrderedCollector.hpp>

#include <mrock/utility/OutputConvenience.hpp>

#include <algorithm>
#include <cmath>
#include <cassert>
#include <string>
#include <iostream>

namespace NickelCUT::commute::Hamiltonian {
using namespace mrock::symbolic_operators;

const std::string outer_K_loop = "#pragma omp parallel for\nfor (int K_pos=0; K_pos < N; ++K_pos) {\nmomentum_iterator<L> K(K_pos);\n";

const std::string independent_interaction_loop =
    "#pragma omp parallel for\n"
    "for (int independent_pos = 0; independent_pos < InteractionDataFrame::independent_size(); ++independent_pos) {\n"
    "const std::size_t tuple_index = InteractionDataFrame::independent_indices()[independent_pos];\n"
    "momentum_iterator<L> K(static_cast<int>(tuple_index / (N * N)));\n"
    "momentum_iterator<L> P(static_cast<int>((tuple_index / N) % N));\n"
    "momentum_iterator<L> Q(static_cast<int>(tuple_index % N));\n";

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

std::string access_expectation_value(const WickOperator& symbol) {
    std::string code;
    if (symbol.type == OperatorType::Number) {
        code += "occupation_numbers";
    }
    else if (symbol.type == OperatorType::HoleNumber) {
        code += "hole_numbers";
    }
    else {
        throw std::invalid_argument("Got an unexpected expectation value.");
    }
    code += "[";
    code += momentum_to_code(symbol.momentum);
    code += "]";
    return code;
}

std::string access_coefficient(const Coefficient& coeff) {
    std::string code = "";
    if (coeff.name == "\\alpha") {
        code += "alpha_sign_cache(";
        code += momentum_to_code(coeff.momenta[0]) + ", " + momentum_to_code(coeff.momenta[1]) + ", " + momentum_to_code(coeff.momenta[2]);
        code += ")\n\t* ";
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
    if (coeff.name == "\\varepsilon") {
        code += "current.dispersion[";
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
    const std::string accessor = "dHdl.dispersion[K] ";

    std::string code = outer_K_loop;
    code += momentum_for_loop("P");
    code += "double Q_sum_value{};\n";
    code += momentum_for_loop("Q");
    
    const IntFractional smallest_factor = get_smallest_factor(bilinears);

    for (const auto& term : bilinears) {
        assert(term.operators[0].momentum == Momentum('P'));
        assert(term.coefficients[0].name == "U");
        assert(term.coefficients[1].name == "\\alpha");

        code += "Q_sum_value ";

        const auto mult = term.multiplicity / smallest_factor;
        code += mult > 0 ? "+= " : "-= ";
        if (mult != 1 && mult != -1) {
            code += std::to_string(std::abs(static_cast<double>(mult)));
            code += " * ";
        }
        
        code += access_coefficient(term.coefficients[0]);
        code += "\n\t* ";
        code += access_coefficient(term.coefficients[1]);
        code += "\n\t";
        // The first expectation value is <n_P> which we drag out of the Q-sum.
        for (auto it = term.operators.begin() + 1; it != term.operators.end(); ++it) {
            code += " * ";
            code += access_expectation_value(*it);
        }
        code += ";\n";
    }

    code += "} // Q-loop\n";
    code += accessor + "+= " + std::to_string(static_cast<double>(smallest_factor)) 
        + " * Q_sum_value * occupation_numbers[P];\n";
    code += "} // P-loop\n";
    code += "} // K-loop\n";
    return code;
}

std::string generate_quartic(const experimental::WickOrderedCollector& quartics, bool parallel) {
    const std::string accessor = std::string("dHdl.") 
        + (parallel ? std::string("interactions_same_spin") : std::string("interactions_differing_spin"))
        + std::string("(K, P, Q) ");

    std::string code = independent_interaction_loop;
    if (parallel) {
        code += "if (K==P) continue; // Pauli principle\n";
    }
    if (parallel) {
        code += "if (K+Q==P-Q) continue; // Pauli principle\n";
    }

    for (auto it = quartics.begin(); it != quartics.end(); ++it) {
        if (!(it->sums.momenta.empty())) continue;
        code += accessor + (it->multiplicity > 0 ? "+= " : "-= ");
        if (it->multiplicity != 1 && it->multiplicity != -1) {
            code += std::to_string(std::abs(static_cast<double>(it->multiplicity)));
            code += " * ";
        }

        code += access_coefficient(it->coefficients[0]);
        code += "\n\t* ";
        code += access_coefficient(it->coefficients[1]);

        assert(it->operators.empty());
        code += ";\n";
    }

    code += momentum_for_loop("R");
    code += "double nR_hQR_value{};\n";

    const auto is_nR_hQR = [](const std::vector<WickOperator>& expecs) {
        if (expecs.size() != 2U) return false;
        if (expecs[0].type != OperatorType::Number || expecs[0].momentum != Momentum('R')) return false;
        return (expecs[1].type == OperatorType::HoleNumber && expecs[1].momentum == Momentum("Q+R"));
    };


    const IntFractional smallest_factor = get_smallest_factor(quartics);

    for (auto& term : quartics) {
        if (term.sums.momenta.empty()) continue;

        code += is_nR_hQR(term.operators) ? "nR_hQR_value " : accessor;

        const auto mult = term.multiplicity / smallest_factor;
        code += mult > 0 ? "+= " : "-= ";
        if (mult != 1 && mult != -1) {
            code += std::to_string(std::abs(static_cast<double>(mult)));
            code += " * ";
        }

        code += access_coefficient(term.coefficients[0]);
        code += "\n\t* ";
        code += access_coefficient(term.coefficients[1]);

        if (!is_nR_hQR(term.operators)) {
            code += "\n\t";
            for (auto it = term.operators.begin(); it != term.operators.end(); ++it) {
                code += " * ";
                code += access_expectation_value(*it);
            }
        }
        code += ";\n";
    }

    code += accessor + "+= nR_hQR_value * occupation_numbers[R] * hole_numbers[Q+R];\n";
    code += "} // R-loop\n";
    code += accessor + "*= " + std::to_string(static_cast<double>(smallest_factor)) + ";\n";

    code += "}\n";

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
        + independent_interaction_loop +
        "#ifdef USE_LW_GENERATOR\n"
        "alpha_sign_cache(K,P,Q) = current.dispersion[K] + current.dispersion[P] - current.dispersion[P-Q] - current.dispersion[K+Q];\n"
        "#else\n"
        "alpha_sign_cache(K,P,Q) = sign(current.dispersion[K] + current.dispersion[P] - current.dispersion[P-Q] - current.dispersion[K+Q]);\n"
        "#endif\n"
        "}\n\n";

    const std::string file_footer = 
        "dHdl.interactions_same_spin.antisymmetrize();\n"
        "dHdl.interactions_same_spin.symmetrize();\n"
        "dHdl.interactions_differing_spin.symmetrize();\n"
        "dHdl.interactions_same_spin.clear_noise();\n"
        "dHdl.interactions_differing_spin.clear_noise();\n"
        "}\n} // namespace NickelCUT::flow";

    const std::string bilinear_code = generate_bilinear(flow_coeffs[0]) 
        + "\n//--------------------------------------------------------------//\n";
    const std::string u_anti_code = generate_quartic(flow_coeffs[1], false)
        + "\n//--------------------------------------------------------------//\n";
    const std::string u_para_code = generate_quartic(flow_coeffs[2], true);

    mrock::utility::save_string_raw(file_header + bilinear_code + u_anti_code + u_para_code + file_footer, "sources/flow/FlowEquation.cpp");
}

}