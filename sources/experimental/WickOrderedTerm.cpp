#include "WickOrderedTerm.hpp"

#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/WickOperator.hpp>
#include <mrock/symbolic_operators/WickOperatorTemplate.hpp>
#include <mrock/symbolic_operators/detail/container_helper.hpp>

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <variant>

namespace mrock::symbolic_operators::experimental {

WickOrderedTerm::WickOrderedTerm(const Term& base) 
    : AbstractTerm<WickOperator>(base.multiplicity,
                                 base.coefficients,
                                 base.sums,
                                 base.delta_momenta,
                                 base.delta_indizes,
                                 std::vector<WickOperator>()),
      wick_expression()
{ }

static WickOrderedTerm from_wick_term_skip_wick_expression(const WickOrderedTerm& base)
{
    WickOrderedTerm result;
    result.multiplicity = base.multiplicity;
    result.coefficients = base.coefficients;
    result.sums = base.sums;
    result.delta_momenta = base.delta_momenta;
    result.delta_indizes = base.delta_indizes;
    result.operators = base.operators;

    return result;
}

void WickOrderedTerm::include_template_result(const TemplateResult::SingleResult& result) 
{
    this->delta_indizes.insert(this->delta_indizes.begin(), result.index_deltas.begin(), result.index_deltas.end());
    this->operators.push_back(result.op);
    this->multiplicity *= result.factor;
}

std::ostream& operator<<(std::ostream& os, const WickOrderedTerm& term)
{
    if (term.multiplicity > 0) {
        os << "+";
    }
    os << term.multiplicity << " ";
    os << term.sums;
    os << term.coefficients << " ";
    for (const auto& delta : term.delta_momenta) {
        os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
    }
    for (const auto& delta : term.delta_indizes) {
        os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
    }
    if (term.is_identity()) {
        os << " \\hat{1} ";
        return os;
    }
    for (const auto& op : term.operators) {
        os << op << " ";
    }
    os << term.wick_expression;

    return os;
}

}