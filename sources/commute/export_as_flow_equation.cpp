#include "export_as_flow_equation.hpp"

#include <mrock/symbolic_operators/WickTerm.hpp>
#include <mrock/utility/OutputConvenience.hpp>

#include <cassert>
#include <string>
#include <iostream>

namespace NickelCUT::commute {
using namespace mrock::symbolic_operators;

std::string momentum_for_loop(std::string const& it_name) {
    return "for (mom_it " + it_name + " = mom_it::begin(); " 
        + it_name + " != mom_it::end(); ++" + it_name + ") {\n";
}

std::string sign_to_string(int factor) {
    return (factor < 0 ? "-" : "+");
}

std::string momentum_to_code(const Momentum& momentum) {
    if (momentum.is_zero()) {
        return "(L + L*L) / 2";
    }

    std::string code = "(";
    for (const auto& symbol : momentum) {
        assert(symbol.factor == 1 || symbol.factor == -1);
        if (code == "(") {
            if (symbol.factor == -1) code += "-";
        }
        else {
            code += sign_to_string(symbol.factor);
        }
        code += symbol.name._n;
    }
    code += ").get_position()";
    return code;
}

std::string access_coefficient(const Coefficient& coeff) {
    std::string code = "";
    if (coeff.name == "\\alpha") {
        code += "sign(";
        
        code += "current.epsilon_tilde[";
        code += momentum_to_code(coeff.momenta[0]);
        code += "] + ";
        code += "current.epsilon_tilde[";
        code += momentum_to_code(coeff.momenta[1]);
        code += "] - ";

        code += "current.epsilon_tilde[";
        code += momentum_to_code(coeff.momenta[0] + coeff.momenta[2]);
        code += "] - ";
        code += "current.epsilon_tilde[";
        code += momentum_to_code(coeff.momenta[1] - coeff.momenta[2]);
        code += "]";

        code += ") \n\t* ";
    }
    if (coeff.name == "U" || coeff.name == "\\alpha") {
        if (coeff.indizes[0] == Index::Parallel) {
            code += "current.interactions_same_spin";
        }
        else {
            code += "current.interactions_differing_spin";
        }
        for (const auto& momentum : coeff.momenta) {
            code += "[";
            code += momentum_to_code(momentum);
            code += "]";
        }
    }
    if (coeff.name == "\\tilde{\\varepsilon}") {
        code += "current.epsilon_tilde[";
        code += momentum_to_code(coeff.momenta[0]);
        code += "]";
    }

    return code;
}

std::string generate_bilinear(const WickTermCollector& bilinears) 
{
    std::string code = momentum_for_loop("p");
    code += "result.dispersion[p.get_position()] = 0.0;\n";
    code += momentum_for_loop("q");
    code += "double nr_value{};\ndouble one_value{};\n";
    code += momentum_for_loop("r");

    for (const auto& term : bilinears) {
        assert(term.operators[0].momentum == Momentum('q'));
        assert(term.coefficients[0].name == "U");
        assert(term.coefficients[1].name == "\\alpha");

        if (term.operators.size() == 1U) {
            code += "one_value += ";
        }
        else {
            assert(term.operators[1].momentum == Momentum('r'));
            code += "nr_value += ";
        }

        code += access_coefficient(term.coefficients[0]);
        code += "\n\t* ";
        code += access_coefficient(term.coefficients[1]);
        code += ";\n";
    }

    code += "nr_value *= current.occupation_numbers[r.get_position()];\n";
    code += "} // r-loop\n";
    code += "result.dispersion[p.get_position()] += (nr_value + one_value) * current.occupation_numbers[q.get_position()];\n";
    code += "} // q-loop\n";
    code += "} // p-loop\n";
    return code;
}



std::string generate_quartic(const WickTermCollector& quartics, bool parallel) {
    const std::string accessor = std::string("result.") 
        + (parallel ? std::string("interactions_same_spin") : std::string("interactions_differing_spin"))
        + std::string("[p.get_position()][q.get_position()][r.get_position()]");

    std::string code = momentum_for_loop("p");
    code += momentum_for_loop("q");
    code += momentum_for_loop("r");
    code += accessor + " = 0.0;\n";

    for (auto it = quartics.begin(); it != quartics.end() && it->sums.momenta.empty(); ++it) {
        code += accessor + " += ";
        code += access_coefficient(it->coefficients[0]);
        code += "\n\t* ";
        code += access_coefficient(it->coefficients[1]);

        if (!(it->operators.empty())) {
            assert(it->operators.size() == 1U);
            code += "\n\t* current.occupation_numbers[";
            code += momentum_to_code(it->operators[0].momentum);
            code += "]";
        }
        code += ";\n";
    }

    code += "double ns_value{};\ndouble one_value{};\n";
    code += momentum_for_loop("s");
    
    for (auto& term : quartics) {
        if (term.sums.momenta.empty()) continue;

        code += (term.operators.empty() ? "one_value += " : "ns_value += ");
        code += access_coefficient(term.coefficients[0]);
        code += "\n\t* ";
        code += access_coefficient(term.coefficients[1]);

        if (term.operators.size() == 2U) {
            code += "\n\t* current.occupation_numbers[";
            if (term.operators[0].momentum == Momentum('s')) {
                code += momentum_to_code(term.operators[1].momentum);
            }
            else if (term.operators[1].momentum == Momentum('s')) {
                code += momentum_to_code(term.operators[0].momentum);
            }
            else {
                code += "QUARTIC ONE SHOULD BE S!";
                return code;
            }
            code += "]";
        }
        code += ";\n";
    }

    code += accessor + " += one_value + current.occupation_numbers[s.get_position()] * ns_value;\n";

    code += "} // s-loop\n";

    code += "} // r-loop\n";
    code += "} // q-loop\n";
    code += "} // p-loop\n";

    return code;
}



void export_as_flow_equation(const std::array<WickTermCollector, 3> flow_coeffs) 
{
    const std::string file_header = "#include \"momentum_iterator.hpp\"\n#include \"FlowContainer.hpp\"\n\n"
        "namespace NickelCUT::flow {\n\n"
        "void flow_step(FlowContainer<L>& result, const FlowContainer<L>& current) {";

    const std::string file_footer = "}\n} // namespace NickelCUT::flow";

    const std::string bilinear_code = generate_bilinear(flow_coeffs[0]) 
        + "\n//--------------------------------------------------------------//\n";
    const std::string u_anti_code = generate_quartic(flow_coeffs[1], false)
        + "\n//--------------------------------------------------------------//\n";
    const std::string u_para_code = generate_quartic(flow_coeffs[2], true);
    
    mrock::utility::save_string_raw(file_header + bilinear_code + u_anti_code + u_para_code + file_footer, "sources/flow/test.cpp");
}

}