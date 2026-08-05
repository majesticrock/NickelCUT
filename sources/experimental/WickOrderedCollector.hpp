#pragma once

#include "WickOrderedTerm.hpp"

#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/WickOperatorTemplate.hpp>
#include <mrock/symbolic_operators/detail/vector_macro.hpp>

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

std::ostream& operator<<(std::ostream& os, const WickOrderedCollector& terms);

} // namespace mrock::symbolic_operators::experimental