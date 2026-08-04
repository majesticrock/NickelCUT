#pragma once
#include "WickOrderedExpression.hpp"

#include <mrock/symbolic_operators/AbstractTerm.hpp>
#include <mrock/symbolic_operators/Operator.hpp>
#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/WickOperator.hpp>
#include <mrock/symbolic_operators/WickOperatorTemplate.hpp>

#include <vector>

namespace mrock::symbolic_operators::experimental {

class WickOrderedTerm : public AbstractTerm<WickOperator> {
public:
    WickOrderedExpression wick_expression;  ///< The Wick ordered expression associated with the term.

    std::vector<Operator> temporary_operators;  ///< Temporary operators used in the term.

    /**
     * @brief Default constructor for WickOrderedTerm.
     */
    WickOrderedTerm() = default;

    /**
     * @brief Constructs a WickOrderedTerm object from a base Term.
     */
    WickOrderedTerm(const Term& base);

    /**
     * @brief Includes a template result in the term.
     *
     * @param result The TemplateResult::SingleResult to include.
     */
    void include_template_result(const TemplateResult::SingleResult& result);
};

std::vector<WickOrderedTerm> wick_decompose(const std::vector<Term>& terms, 
    const std::vector<WickOperatorTemplate>& templates);

} // namespace mrock::symbolic_operators::experimental