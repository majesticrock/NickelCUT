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
    std::sort(operators.begin(), operators.end());
    std::sort(sums.momenta.begin(), sums.momenta.end());
    std::sort(sums.spins.begin(), sums.spins.end());
}

bool WickOrderedTerm::is_identity() const noexcept {
    return wick_expression.empty();
}

bool WickOrderedTerm::is_bilinear() const noexcept {
    return wick_expression.size() == 2U;
}

bool WickOrderedTerm::is_quartic() const noexcept {
    return wick_expression.size() == 4U;
}

WickOrderedTerm& WickOrderedTerm::hermitian_conjugate_inplace() noexcept
{
    wick_expression.hermitian_conjugate_inplace();
    return *this;
}

WickOrderedTerm WickOrderedTerm::hermitian_conjugate() const noexcept
{
    WickOrderedTerm copy(*this);
    copy.hermitian_conjugate_inplace();
    return copy;
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