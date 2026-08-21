#pragma once

#include "WickOrderedTerm.hpp"

#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/WickOperatorTemplate.hpp>
#include <mrock/symbolic_operators/WickSymmetry.hpp>
#include <mrock/symbolic_operators/detail/vector_macro.hpp>

#include <memory>
#include <ostream>
#include <vector>


namespace mrock::symbolic_operators::experimental {

/**
 * @brief A collector for managing WickOrderedTerm objects.
 */
struct WickOrderedCollector {
    std::vector<WickOrderedTerm> terms;  ///< The collected \c WickTerm objects

    /**
     * @brief Default constructor for WickOrderedCollector.
     */
    WickOrderedCollector() = default;

    /**
     * @brief Combines duplicate terms
     */
    void combine_duplicates();

    /**
     * @brief Takes the Hermitian conjugate of the current expression in place
     * 
     * @return A reference to \c *this after the Hermitian conjugation 
     */
    WickOrderedCollector& hermitian_conjugate_inplace() noexcept;

    /**
     * @brief Takes the Hermitian conjugate of the current expression in a copy
     * 
     * @return A copy of \c *this after the Hermitian conjugation 
     */
    WickOrderedCollector hermitian_conjugate() const noexcept;

    WickOrderedCollector& operator+=(const WickOrderedCollector& other);

    WickOrderedCollector& operator-=(const WickOrderedCollector& other);

    MROCK_VECTOR_WRAPPER_FILL_MEMBERS(WickOrderedTerm, terms);
};

/**
 * @brief Decomposes a Term into WickOrderedTerm objects using the provided WickOperatorTemplates.
 * 
 * @param term The Term to decompose.
 * @param templates The vector of WickOperatorTemplate objects to use for decomposition.
 * @return WickOrderedCollector The collector containing the resulting WickOrderedTerm objects.
 */
WickOrderedCollector wick_decompose(const Term& term, const std::vector<WickOperatorTemplate>& templates);

WickOrderedCollector wick_decompose(const std::vector<Term>& terms, const std::vector<WickOperatorTemplate>& templates);

void clean_wick_ordered_terms(WickOrderedCollector& terms,
    const std::vector<std::unique_ptr<WickSymmetry>>& symmetries = std::vector<std::unique_ptr<WickSymmetry>>{});

std::ostream& operator<<(std::ostream& os, const WickOrderedCollector& terms);

WickOrderedCollector operator+(WickOrderedCollector lhs, const WickOrderedCollector& rhs);

WickOrderedCollector operator-(WickOrderedCollector lhs, const WickOrderedCollector& rhs);

} // namespace mrock::symbolic_operators::experimental