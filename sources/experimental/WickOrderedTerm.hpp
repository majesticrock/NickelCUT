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
};

std::ostream& operator<<(std::ostream& os, const WickOrderedTerm& term);

} // namespace mrock::symbolic_operators::experimental