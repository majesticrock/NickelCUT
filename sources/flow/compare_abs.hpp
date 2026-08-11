#pragma once

#include <cmath>

namespace NickelCUT::flow
{
struct LessThanAbs {
    constexpr bool operator()(const double lhs, const double rhs) {
        return std::abs(lhs) < std::abs(rhs);
    }  
};
} // namespace NickelCUT::flow
