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
void WickOrderedTerm::for_each_momentum_except_deltas(const std::function<void(Momentum&)>& f) 
{
    AbstractTerm<WickOperator>::for_each_momentum_except_deltas(f);
    for (auto& op : wick_expression) {
        f(op.momentum);
    }
}

void WickOrderedTerm::for_each_index_except_deltas(const std::function<void(IndexWrapper&)>& f) 
{
    AbstractTerm<WickOperator>::for_each_index_except_deltas(f);
    for (auto& op : wick_expression) {
        f(op.indizes);
    }
}

WickOrderedTerm::WickOrderedTerm(const Term& base) 
    : AbstractTerm<WickOperator>(base.multiplicity,
                                 base.coefficients,
                                 base.sums,
                                 base.delta_momenta,
                                 base.delta_indizes,
                                 std::vector<WickOperator>()),
      wick_expression()
{ }

WickOrderedTerm WickOrderedTerm::from_wick_term_skip_wick_expression(const WickOrderedTerm& base)
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

bool WickOrderedTerm::resolve_deltas()
{
    if (!resolve_momentum_deltas())
        return false;
    if (!resolve_index_deltas())
        return false;

    return true;
}

void WickOrderedTerm::discard_zero_momenta() 
{
    AbstractTerm<WickOperator>::discard_zero_momenta();
    for (auto& op : wick_expression) {
        op.momentum.remove_zeros();
    }
}

void WickOrderedTerm::sort()
{
    for (auto& delta : delta_momenta) {
        if (delta.first.momentum_list.size() == 1 && delta.second.momentum_list.size() == 1) {
            // This comparison is well defined because we save the momentum as char i.e. byte
            // which is easily comparable
            if (delta.first.momentum_list[0].name < delta.second.momentum_list[0].name) {
                std::swap(delta.first, delta.second);
                if (delta.first.momentum_list[0].factor < 0) {
                    delta.first.flip_momentum();
                    delta.second.flip_momentum();
                }
                if (delta.first.add_Q) {
                    delta.first.add_Q = false;
                    delta.second.add_Q = !(delta.second.add_Q);
                }
            }
            for (auto& op : operators) {
                op.momentum.replace_occurances(delta.first.momentum_list[0].name, delta.second);
            }
            for (auto& coeff : coefficients) {
                coeff.momenta.replace_occurances(delta.first.momentum_list[0].name, delta.second);
            }
        }
    }

    for (auto& op : operators) {
        if (op.type == OperatorType::CDW && op.momentum.add_Q) {
            op.momentum.add_Q = false;
            op.is_daggered = !(op.is_daggered);
        }
    }

    for (std::size_t i = 0U; i < operators.size(); ++i) {
        for (std::size_t j = i + 1U; j < operators.size(); ++j) {
            if (operators[i].type > operators[j].type) {
                std::swap(operators[i], operators[j]);
            } else if (operators[i].type == operators[j].type) {
                if (momentum_order(operators[i].momentum, operators[j].momentum)) {
                    std::swap(operators[i], operators[j]);
                }
            }
        }
    }

    for (auto& coeff : coefficients) {
        for (auto& momentum : coeff.momenta) {
            momentum.sort();

            if (coeff.inversion_symmetry && !momentum.momentum_list.empty()) {
                if (momentum.momentum_list[0].factor < 0) {
                    momentum.flip_momentum();
                }
            }
            if (coeff.Q_changes_sign && momentum.add_Q) {
                momentum.add_Q = false;
                this->multiplicity *= -1;
            }
        }
    }

    for (auto& coeff : coefficients) {
        if (coeff.momenta.size() == 3U) {
            Momentum* first_momentum = &coeff.momenta.front();
            if (first_momentum != nullptr && first_momentum->empty()) {
                if (coeff.momenta.size() > 1U)
                    first_momentum = &coeff.momenta[1];
            }
            if ((first_momentum != nullptr) && (!first_momentum->empty())) {
                if (!first_momentum->first_momentum_is('k')) {
                    coeff.use_symmetric_interaction_exchange();
                }
                if (sums.momenta.empty()) {
                    if (coeff.momenta.back().first_momentum_is_negative()) {
                        coeff.use_symmetric_interaction_inversion();
                    }
                } else if ((!sums.momenta.is_summed_over(first_momentum->front().name)) &&
                           first_momentum->first_momentum_is_negative()) {
                    coeff.use_symmetric_interaction_inversion();
                }
            }
        }

        for (auto& momentum : coeff.momenta) {
            if (momentum.empty())
                continue;
            if ((!sums.momenta.is_summed_over(momentum.front().name)) && momentum.front().factor < 0) {
                if (coeff.inversion_symmetry)
                    momentum.flip_momentum();
            }

            for (const auto& sum : sums.momenta) {
                int idx = momentum.is_used_at(sum);
                if (idx < 0)
                    continue;

                if (momentum.momentum_list[idx].factor < 0) {
                    invert_momentum_sum(sum);
                }
            }
        }
    }
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
    for (const auto& op : term.operators) {
        os << op << " ";
    }
    os << term.wick_expression;

    return os;
}

bool operator==(const WickOrderedTerm& lhs, const WickOrderedTerm& rhs) {
    if (lhs.coefficients != rhs.coefficients)
        return false;
    if (lhs.sums != rhs.sums)
        return false;
    if (lhs.delta_indizes != rhs.delta_indizes)
        return false;
    if (lhs.delta_momenta != rhs.delta_momenta)
        return false;
    if (lhs.operators != rhs.operators)
        return false;
    if (lhs.wick_expression != rhs.wick_expression)
        return false;
    return true;
}

bool operator!=(const WickOrderedTerm& lhs, const WickOrderedTerm& rhs) {
    return !(lhs == rhs);
}
}