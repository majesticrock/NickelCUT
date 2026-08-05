#include "WickOrderedExpression.hpp"
#include <ostream>

namespace mrock::symbolic_operators::experimental {

std::ostream& operator<<(std::ostream& os, const WickOrderedExpression& expr) {
    if (expr.empty()) {
        os << ": \\hat{1} :";
    }
    else {
        os << ":" << expr.operators << ":";
    }
    return os;
}

bool operator==(const WickOrderedExpression& lhs, const WickOrderedExpression& rhs) {
    return lhs.operators == rhs.operators;
}

bool operator!=(const WickOrderedExpression& lhs, const WickOrderedExpression& rhs) {
    return !(lhs == rhs);
}
}