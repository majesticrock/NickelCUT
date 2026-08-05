#include "WickOrderedCollector.hpp"

#include <utility>

namespace mrock::symbolic_operators::experimental {

// Computes n!!
template <typename IntegerType, typename RealType = std::size_t>
constexpr RealType double_factorial(const IntegerType n) {
    RealType result{1};
    for (std::make_signed_t<IntegerType> i = n; i > 0; i -= 2) {
        result *= i;
    }
    return result;
}

void wick_decompose_internal(const std::vector<Operator>& operators, 
    const std::vector<WickOperatorTemplate>& templates,
    WickOrderedTerm&& attach_to,
    WickOrderedCollector& result)
{
    result.emplace_back(std::move(attach_to));
    result.back().wick_expression.operators = operators;

    if (operators.size() < 2U) {
        return;
    }

    for (std::size_t i = 1U; i < operators.size(); ++i) {
        std::vector<Operator> copy_operators = operators;
        copy_operators.erase(copy_operators.begin() + i);
        copy_operators.erase(copy_operators.begin());

        for (const auto& operator_template : templates) {
            auto template_result = operator_template.create_from_operators(operators[0], operators[i]);

            if (template_result) {
                for (const auto& single_result : template_result.results) {
                    if (single_result.contains_impossible_delta()) {
                        continue;
                    }

                    WickOrderedTerm new_term = WickOrderedTerm::from_wick_term_skip_wick_expression(result.back());
                    new_term.wick_expression.operators.clear();
                    if (!(i & 1)) {
                        new_term.multiplicity *= -1;
                    }

                    new_term.include_template_result(single_result);
                    wick_decompose_internal(copy_operators, templates, std::move(new_term), result);
                }                
            }
        }
    }
}


WickOrderedCollector wick_decompose(const Term& term, 
    const std::vector<WickOperatorTemplate>& templates)
{
    const std::size_t estimated_size = double_factorial(term.get_operators().size()) * templates.size();
    WickOrderedCollector result;
    result.reserve(estimated_size);

    WickOrderedTerm initial(term);

    wick_decompose_internal(term.get_operators(), templates, std::move(initial), result);

    return result;
}

std::ostream& operator<<(std::ostream& os, const WickOrderedCollector& terms) {
    for (WickOrderedCollector::const_iterator it = terms.begin(); it != terms.end(); ++it) {
        os << "\t&" << *it;
        if (it != terms.end() - 1) {
            os << " \\\\";
        }
        os << "\n";
    }
    return os;
}
}