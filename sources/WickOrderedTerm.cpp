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

// Computes n!!
template <typename IntegerType, typename RealType = std::size_t>
constexpr RealType double_factorial(const IntegerType n) {
    RealType result{1};
    for (std::make_signed_t<IntegerType> i = n; i > 0; i -= 2) {
        result *= i;
    }
    return result;
}

WickOrderedTerm::WickOrderedTerm(const Term& base)
    : AbstractTerm<WickOperator>(base.multiplicity,
                                 base.coefficients,
                                 base.sums,
                                 base.delta_momenta,
                                 base.delta_indizes,
                                 std::vector<WickOperator>()),
        wick_expression(), temporary_operators() 
{}

void WickOrderedTerm::include_template_result(const TemplateResult::SingleResult& result) 
{
    // TODO: Not yet done; Missing logic for wick ordered expressions
    this->delta_indizes.insert(this->delta_indizes.begin(), result.index_deltas.begin(), result.index_deltas.end());
    this->operators.push_back(result.op);
    this->multiplicity *= result.factor;
}

void wick_processor(const std::vector<Operator>& remaining,
                    std::vector<WickOrderedTerm>& reciever_list,
                    std::variant<WickOrderedTerm, Term> buffer) 
{
    if (remaining.empty()) {
        reciever_list.push_back(std::get<WickOrderedTerm>(buffer));
        return;
    }
    for (std::size_t i = 1U; i < remaining.size(); ++i) {
        if (std::holds_alternative<Term>(buffer)) {
            WickOrderedTerm temp(std::get<Term>(buffer));
            buffer = temp;
        }
        if (!(i & 1)) {
            std::get<WickOrderedTerm>(buffer).multiplicity *= -1;
        }
        std::get<WickOrderedTerm>(buffer).temporary_operators.reserve(
            std::get<WickOrderedTerm>(buffer).temporary_operators.size() + 2);
        std::get<WickOrderedTerm>(buffer).temporary_operators.push_back(remaining[0]);
        std::get<WickOrderedTerm>(buffer).temporary_operators.push_back(remaining[i]);

        std::vector<Operator> copy_operators = remaining;
        copy_operators.erase(copy_operators.begin() + i);
        copy_operators.erase(copy_operators.begin());
        wick_processor(copy_operators, reciever_list, buffer);

        // delete last two elements, as they are to be updated in the next iteration
        std::get<WickOrderedTerm>(buffer).temporary_operators.pop_back();
        std::get<WickOrderedTerm>(buffer).temporary_operators.pop_back();
        if (!(i & 1)) {
            std::get<WickOrderedTerm>(buffer).multiplicity *= -1;
        }
    }
}


std::vector<WickOrderedTerm> identify_wick_operators(const WickOrderedTerm& source,
                                          const std::vector<WickOperatorTemplate>& operator_templates) {
    std::vector<WickOrderedTerm> ret;
    ret.push_back(source);
    ret.back().temporary_operators.clear();

    for (std::size_t i = 0U; i < source.temporary_operators.size(); i += 2U) {
        std::vector<TemplateResult> template_results;
        for (const auto& operator_template : operator_templates) {
            auto template_result = operator_template.create_from_operators(source.temporary_operators[i],
                                                                           source.temporary_operators[i + 1U]);
            if (template_result) {
                template_results.push_back(std::move(template_result));
            }
        }

        const std::size_t current_size = ret.size();
        const std::size_t number_additional_elements =
            std::accumulate(template_results.begin(), template_results.end(), std::size_t{},
                            [](std::size_t current, const TemplateResult& tr) { return current + tr.results.size(); });
        if (number_additional_elements > 1U) {
            duplicate_n_inplace(ret, number_additional_elements - 1U);
        }

        std::size_t template_result_it{};
        std::size_t old_it{};
        for (const auto& tr : template_results) {
            old_it = template_result_it;
            for (const auto& tr_result : tr.results) {
                for (auto it = ret.begin(); it != ret.begin() + current_size; ++it) {
                    (it + template_result_it * current_size)->include_template_result(tr_result);
                }
                ++template_result_it;
            }
            std::for_each(ret.begin() + old_it * current_size, ret.begin() + current_size * template_result_it,
                          [&tr](WickOrderedTerm& ret_element) {
                              if (!tr.momentum_delta.isOne())
                                  ret_element.delta_momenta.push_back(tr.momentum_delta);
                          });
        }
    }

    return ret;
}

std::vector<WickOrderedTerm> wick_decompose(const std::vector<Term>& terms, 
    const std::vector<WickOperatorTemplate>& templates)
{
    const std::size_t estimated_size = std::accumulate(
        terms.begin(), terms.end(), std::size_t{},
        [](std::size_t current, const Term& term) { return current + double_factorial(term.get_operators().size()); });

    std::vector<WickOrderedTerm> prepared_wick;
    prepared_wick.reserve(estimated_size);

    for (const auto& term : terms) {
        if (term.is_identity()) {
            prepared_wick.push_back(WickOrderedTerm(term));
        } else {
            wick_processor(term.get_operators(), prepared_wick, term);
        }
    }

    std::vector<WickOrderedTerm> result;
    result.reserve(prepared_wick.size() * templates.size());

    for (auto& w_term : prepared_wick) {
        append_if(result, identify_wick_operators(w_term, templates), [](const WickOrderedTerm& wick) {
            return !(is_always_zero(wick.delta_indizes) || is_always_zero(wick.delta_momenta));
        });
    }
    return result;
}

}