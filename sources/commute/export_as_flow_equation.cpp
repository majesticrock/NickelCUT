#include "export_as_flow_equation.hpp"

#include <mrock/symbolic_operators/WickTerm.hpp>

#include <cassert>
#include <string>


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
    std::string code = "";
    for (const auto& symbol : momentum) {
        assert(symbol.factor == 1 || symbol.factor == -1);
        if (code == "" && symbol.factor == -1) {
            code += "-" + symbol.name._n;
        }
        else {
            code += sign_to_string(symbol.factor) + symbol.name._n;
        }
    }
    return code;
}

std::string access_coefficient(const Coefficient& coeff) {
    if (coeff.name == "\\alpha") {

    }
    if (coeff.name == "U" || coeff.name == "\\alpha") {

    }
    if (coeff.name == "\\tilde{\\varepsilon}") {

    }

    return "INVALID COEFFICIENT";
}

std::string generate_bilinear(const WickTermCollector& bilinears) 
{
    std::string code = momentum_for_loop("p");
    code += "double value{};\n";
    code += momentum_for_loop("q");
    code += "double nr_value{};\ndouble one_value{};\n";
    code += momentum_for_loop("r");

    for (const auto& term : bilinears) {
        assert(term.operators[0].momentum == Momentum('q'));
        if (term.operators.size() == 1U) {
            assert(term.coefficients[0].name == "U");
            assert(term.coefficients[1].name == "\\alpha");

            if(term.sums.spins.empty()) {
                assert(term.coefficients[0].indizes[0] == term.coefficients[0].indizes[1]);
                assert(term.coefficients[1].indizes[0] == term.coefficients[1].indizes[1]);
                code += "one_value += current.interactions_same_spin[(";
                code += ")]";
            }
            else {
                assert(term.coefficients[0].indizes[0] != term.coefficients[0].indizes[1]);
                assert(term.coefficients[1].indizes[0] != term.coefficients[1].indizes[1]);
                code += "one_value += ";
            }
            
        }
        else {
            assert(term.operators[1].momentum == Momentum('r'));
        }
    }

    code += "}\n"; // r-loop
    code += "}\n"; // q-loop
    code += "value *= current.occupation_numbers[q.get_position()];\n";
    code += "}"; // p-loop
    return code;
}

void NickelCUT::commute::export_as_flow_equation(const WickTermCollector& bilinears,
                                                 const WickTermCollector& quartics) 
{

}

}