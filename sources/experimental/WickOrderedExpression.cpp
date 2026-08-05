#include "WickOrderedExpression.hpp"
#include <ostream>

namespace mrock::symbolic_operators::experimental {

std::ostream& operator<<(std::ostream& os, const WickOrderedExpression& expr) {
    os << ":" << expr.operators << ":";
    return os;
}

}