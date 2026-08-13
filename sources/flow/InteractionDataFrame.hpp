#pragma once

#include "compare_abs.hpp"
#include "momentum_iterator.hpp"

#include <cassert>
#include <utility>
#include <algorithm>
#include <array>

namespace NickelCUT::flow {

template<int _N>
class InteractionDataFrame {
    static constexpr int total_size = _N * _N * _N;
    std::array<double, total_size> _data;

public:
    constexpr InteractionDataFrame() {
        _data.fill(double{});
    }
    constexpr InteractionDataFrame(double value) {
        _data.fill(value);
    }
    constexpr InteractionDataFrame(const std::array<std::array<std::array<double, _N>, _N>, _N>& data)
        : _data(data) {};
    constexpr InteractionDataFrame(std::array<std::array<std::array<double, _N>, _N>, _N>&& data)
        : _data(std::move(data)) {};

    constexpr double& operator()(std::size_t x,
                                 std::size_t y,
                                 std::size_t z) noexcept {
        assert(x < _N);
        assert(y < _N);
        assert(z < _N);
        return _data[z + _N * (y + _N * x)];
    }

    constexpr const double& operator()(std::size_t x,
                                       std::size_t y,
                                       std::size_t z) const noexcept {
        assert(x < _N);
        assert(y < _N);
        assert(z < _N);
        return _data[z + _N * (y + _N * x)];
    }

    constexpr void symmetrize() noexcept {
        for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
        for (momentum_iterator<L> q = momentum_iterator<L>::begin(); q != momentum_iterator<L>::end(); ++q) {
        for (momentum_iterator<L> r = momentum_iterator<L>::begin(); r != momentum_iterator<L>::end(); ++r) {
            (*this)(p.get_position(), q.get_position(), r.get_position()) = 0.5 * (
                (*this)(p.get_position(), q.get_position(), r.get_position()) 
                + (*this)(q.get_position(), p.get_position(), (-r).get_position())
            );
            (*this)(q.get_position(), p.get_position(), (-r).get_position()) = (*this)(p.get_position(), q.get_position(), r.get_position());
        }}}
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
        for (std::size_t i=0U; i < total_size; ++i) {
            _data[i] += other._data[i];
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator-=(const InteractionDataFrame& other) noexcept {
        for (std::size_t i=0U; i < total_size; ++i) {
            _data[i] -= other._data[i];
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator*=(const InteractionDataFrame& other) noexcept {
        for (std::size_t i=0U; i < total_size; ++i) {
            _data[i] *= other._data[i];
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator/=(const InteractionDataFrame& other) noexcept {
        for (std::size_t i=0U; i < total_size; ++i) {
            _data[i] /= other._data[i];
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator*=(const double other) noexcept {
        for (std::size_t i=0U; i < total_size; ++i) {
            _data[i] *= other;
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator/=(const double other) noexcept {
        for (std::size_t i=0U; i < total_size; ++i) {
            _data[i] /= other;
        }
        return *this;
    }

    constexpr std::array<double, total_size>& get_data() noexcept {
        return _data;
    }

    constexpr const std::array<double, total_size>& get_data() const noexcept {
        return _data;
    }

    constexpr std::array<std::array<std::array<double, _N>, _N>, _N> as_3D_array() const noexcept {
        std::array<std::array<std::array<double, _N>, _N>, _N> result;

        for (std::size_t x=0; x < _N; ++x) {
            for (std::size_t y=0; y < _N; ++y) {
                for (std::size_t z=0; z < _N; ++z) {
                    result[x][y][z] = (*this)(x, y, z);
                }
            }
        }
        return result;
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