#pragma once
#include "WickOrderedExpression.hpp"

#include <mrock/symbolic_operators/AbstractTerm.hpp>
#include <mrock/symbolic_operators/Operator.hpp>
#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/WickOperator.hpp>
#include <mrock/symbolic_operators/WickOperatorTemplate.hpp>

#include <ostream>
#include <vector>

namespace mrock::symbolic_operators::experimental {

class WickOrderedTerm : public AbstractTerm<WickOperator> {
protected:
    void for_each_momentum_except_deltas(const std::function<void(Momentum&)>& f) override;

    void for_each_index_except_deltas(const std::function<void(IndexWrapper&)>& f) override;

public:
    WickOrderedExpression wick_expression;  ///< The Wick ordered expression associated with the term.
    
    /**
     * @brief Default constructor for WickOrderedTerm.
     */
    WickOrderedTerm() = default;

    /**
     * @brief Constructs a WickOrderedTerm object using the coefficients, sums, etc provided from a base Term
     * 
     * @param base The base Term object.
     */
    WickOrderedTerm(const Term& base);

    /**
     * @brief Constructs a WickOrderedTerm object using the coefficients, sums, etc provided from a base WickTerm, while ignoring the Wick expression of the base term.
     * 
     * @param base The base WickTerm object.
     */
    static WickOrderedTerm from_wick_term_skip_wick_expression(const WickOrderedTerm& base);

    /**
     * @brief Includes a template result in the term.
     *
     * @param result The TemplateResult::SingleResult to include.
     */
    void include_template_result(const TemplateResult::SingleResult& result);

    /**
     * @brief Resolves the deltas in the term.
     */
    bool resolve_deltas();

    /**
     * @brief Discards momenta that are zero from the term.
     */
    void discard_zero_momenta();

    /**
     * @brief Sorts the operators in the term.
     */
    void sort();
};

/**
 * @brief Overloads the output stream operator for WickOrderedTerm.
 * 
 * @param os The output stream.
 * @param term The WickOrderedTerm to output.
 * @return The output stream.
 */
std::ostream& operator<<(std::ostream& os, const WickOrderedTerm& term);

/**
 * @brief Compares two WickOrderedTerm objects for equality.
 * 
 * @param lhs The left-hand side WickOrderedTerm.
 * @param rhs The right-hand side WickOrderedTerm.
 * @return true if the two WickOrderedTerm objects are equal, false otherwise.
 */
bool operator==(const WickOrderedTerm& lhs, const WickOrderedTerm& rhs);

/**
 * @brief Compares two WickOrderedTerm objects for inequality.
 * 
 * @param lhs The left-hand side WickOrderedTerm.
 * @param rhs The right-hand side WickOrderedTerm.
 * @return true if the two WickOrderedTerm objects are not equal, false otherwise.
 */
bool operator!=(const WickOrderedTerm& lhs, const WickOrderedTerm& rhs);

} // namespace mrock::symbolic_operators::experimental