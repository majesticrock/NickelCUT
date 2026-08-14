#pragma once

#include <limits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <bit>

namespace NickelCUT {
/* This function abuses the structure of our desired precision:
 * The mantissa is empty, i.e., we can solely rely on the exponent.
 * If the exponent of <number> is >= 0b01111010011, |number| >= precision, i.e., not 0
 */
constexpr bool is_zero(double number) {
    static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 floating point not verified!");
    // 5.684341886080802e-14   <->   0 | 01111010011 | 0000000000000000000000000000000000000000000000000000
    const uint64_t bits = std::bit_cast<uint64_t>(number);
    return static_cast<uint16_t>((bits >> 52) & 0x7FF) < 0b01111010011;
}

constexpr double sign(double x) {
    if (is_zero(x)) return 0.0;
    return (x < 0.0 ? -1. : 1.);
}

constexpr bool float_equal(double l, double r) {
    return is_zero(std::abs(l - r));
}

struct LessThanAbs {
    constexpr bool operator()(const double lhs, const double rhs) {
        return std::abs(lhs) < std::abs(rhs);
    }  
};

}