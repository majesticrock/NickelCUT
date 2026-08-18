#pragma once

#include <limits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <bit>

namespace NickelCUT {

constexpr std::uint16_t EXPONENT_OF_THRESHOLD = 0b01111011000u;

/* This function abuses the structure of our desired precision:
 * The mantissa is empty, i.e., we can solely rely on the exponent.
 * If the exponent of <number> is >= EXPONENT_OF_THRESHOLD, |number| >= precision, i.e., not 0
 */
constexpr bool is_zero(double number) noexcept {
    static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 floating point not verified!");
    // 1.8189894035458565E-12   <->   0 | 01111011000 | 0000000000000000000000000000000000000000000000000000
    const uint64_t bits = std::bit_cast<uint64_t>(number);
    return static_cast<uint16_t>((bits >> 52) & 0x7FF) < EXPONENT_OF_THRESHOLD;
}
constexpr double exponent_to_double(std::uint16_t exponent) noexcept {
    return std::bit_cast<double>(
        static_cast<std::uint64_t>(exponent) << 52
    );
}
constexpr double PRECISION = exponent_to_double(EXPONENT_OF_THRESHOLD);

constexpr auto threshold =
    exponent_to_double(0b01111010011);

constexpr double sign(double x) noexcept {
    if (is_zero(x)) return 0.0;
    return (x < 0.0 ? -1. : 1.);
}

constexpr bool float_equal(double l, double r) noexcept {
    return is_zero(l - r);
}

struct LessThanAbs {
    constexpr bool operator()(const double lhs, const double rhs) const noexcept {
        return std::abs(lhs) < std::abs(rhs);
    }  
};

struct LessOrAlmostEqual {
    constexpr bool operator()(const double lhs, const double rhs) const noexcept {
        if (lhs < rhs) return true;
        return is_zero(lhs - rhs);
    }
};

struct GreaterOrAlmostEqual {
    constexpr bool operator()(const double lhs, const double rhs) const noexcept {
        if (lhs > rhs) return true;
        return is_zero(lhs - rhs);
    }
};

constexpr bool less_or_almost_equal(const double lhs, const double rhs) noexcept {
    return LessOrAlmostEqual()(lhs, rhs);
}

constexpr bool greater_or_almost_equal(const double lhs, const double rhs) noexcept {
    return GreaterOrAlmostEqual()(lhs, rhs);
}


}