#pragma once

#include "compare_abs.hpp"

#include <utility>
#include <algorithm>
#include <array>

namespace NickelCUT::flow {

template<int _N>
class InteractionDataFrame {
    std::array<double, _N*_N*_N> _data;

public:
    constexpr InteractionDataFrame() {
        _data.fill(double{});
    };
    constexpr InteractionDataFrame(double value) {
        _data.fill(value);
    };
    constexpr InteractionDataFrame(const std::array<std::array<std::array<double, _N>, _N>, _N>& data)
        : _data(data) {};
    constexpr InteractionDataFrame(std::array<std::array<std::array<double, _N>, _N>, _N>&& data)
        : _data(std::move(data)) {};

    constexpr double& operator()(std::size_t x,
                                 std::size_t y,
                                 std::size_t z) noexcept {
        return _data[z + _N * (y + _N * x)];
    }

    constexpr const double& operator()(std::size_t x,
                                       std::size_t y,
                                       std::size_t z) const noexcept {
        return _data[z + _N * (y + _N * x)];
    }

    // Required for boost odeint
    constexpr double abs() const noexcept {  ///< uses the L1 norm
        double val{};
        for(const auto& element : _data) {
            val += std::abs(element);
        }
        return val;
    };
    constexpr double norm_inf() const noexcept {  ///< uses the L_infinity norm
        return *std::max_element(_data.begin(), _data.end(), LessThanAbs());
    };

    constexpr InteractionDataFrame& operator+=(const InteractionDataFrame& other) noexcept {
        for (std::size_t x = 0U; x < _N; ++x) {
            for (std::size_t y = 0U; y < _N; ++y) {
                for (std::size_t z = 0U; z < _N; ++z) {
                    (*this)(x, y, z) += other(x, y, z);
                }
            }
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator-=(const InteractionDataFrame& other) noexcept {
        for (std::size_t x = 0U; x < _N; ++x) {
            for (std::size_t y = 0U; y < _N; ++y) {
                for (std::size_t z = 0U; z < _N; ++z) {
                    (*this)(x, y, z) -= other(x, y, z);
                }
            }
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator*=(const InteractionDataFrame& other) noexcept {
        for (std::size_t x = 0U; x < _N; ++x) {
            for (std::size_t y = 0U; y < _N; ++y) {
                for (std::size_t z = 0U; z < _N; ++z) {
                    (*this)(x, y, z) *= other(x, y, z);
                }
            }
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator/=(const InteractionDataFrame& other) noexcept {
        for (std::size_t x = 0U; x < _N; ++x) {
            for (std::size_t y = 0U; y < _N; ++y) {
                for (std::size_t z = 0U; z < _N; ++z) {
                    (*this)(x, y, z) /= other(x, y, z);
                }
            }
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator*=(const double other) noexcept {
        for (std::size_t x = 0U; x < _N; ++x) {
            for (std::size_t y = 0U; y < _N; ++y) {
                for (std::size_t z = 0U; z < _N; ++z) {
                    (*this)(x, y, z) *= other;
                }
            }
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator/=(const double other) noexcept {
        for (std::size_t x = 0U; x < _N; ++x) {
            for (std::size_t y = 0U; y < _N; ++y) {
                for (std::size_t z = 0U; z < _N; ++z) {
                    (*this)(x, y, z) /= other;
                }
            }
        }
        return *this;
    }

    constexpr std::array<double, _N*_N*_N>& get_data() noexcept {
        return _data;
    }

    constexpr const std::array<double, _N*_N*_N>& get_data() const noexcept {
        return _data;
    }
};

template <int _N>
constexpr InteractionDataFrame<_N> operator+(InteractionDataFrame<_N> lhs, const InteractionDataFrame<_N>& rhs) { return (lhs += rhs); }
template <int _N>
constexpr InteractionDataFrame<_N> operator-(InteractionDataFrame<_N> lhs, const InteractionDataFrame<_N>& rhs) { return (lhs -= rhs); }
template <int _N>
constexpr InteractionDataFrame<_N> operator*(InteractionDataFrame<_N> lhs, const InteractionDataFrame<_N>& rhs) { return (lhs *= rhs); }
template <int _N>
constexpr InteractionDataFrame<_N> operator/(InteractionDataFrame<_N> lhs, const InteractionDataFrame<_N>& rhs) { return (lhs /= rhs); }

template <int _N>
constexpr InteractionDataFrame<_N> operator*(InteractionDataFrame<_N> lhs, const double rhs) { return (lhs *= rhs); }
template <int _N>
constexpr InteractionDataFrame<_N> operator*(const double lhs, InteractionDataFrame<_N> rhs) { return (rhs *= lhs); }
template <int _N>
constexpr InteractionDataFrame<_N> operator/(InteractionDataFrame<_N> lhs, const double rhs) { return (lhs /= rhs); }

} // namespace NickelCut::flow