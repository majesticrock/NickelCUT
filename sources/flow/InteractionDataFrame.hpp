#pragma once

#include "../helper_functions.hpp"
#include "momentum_iterator.hpp"

#include <omp.h>

#include <cassert>
#include <cstddef>
#include <utility>
#include <algorithm>
#include <vector>
#include <iostream>

namespace NickelCUT::flow {

template<int _N>
class InteractionDataFrame {
    static constexpr int total_size = _N * _N * _N;
    std::vector<double> _data;

public:
    template<class Archive>
    void serialize(Archive& ar, [[maybe_unused]] const unsigned int version) {
        ar & _data;
    }

    InteractionDataFrame()
        : _data(total_size, 0.0)
    {}
    InteractionDataFrame(double value)
        : _data(total_size, value)
    {}

    constexpr InteractionDataFrame(const std::vector<double>& data)
        : _data(data) { assert(data.size() == total_size); };
    constexpr InteractionDataFrame(std::vector<double>&& data)
        : _data(std::move(data)) { assert(data.size() == total_size); };

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

    constexpr void fill(double value) noexcept {
        for (auto& element : _data) {
            element = value;
        }
    }

    constexpr double abs_squared() const noexcept {  ///< the square of the L2 norm
        double val{};
        for(const auto& element : _data) {
            val += element*element;
        }
        return val;
    };

    constexpr std::size_t size() const noexcept {
        return _data.size();
    }

    // Required for boost odeint
    constexpr double abs() const noexcept {  ///< the L2 norm
        return std::sqrt(abs_squared());
    };
    constexpr double norm_inf() const noexcept {  ///< the L_infinity norm
        return std::abs(*std::max_element(_data.begin(), _data.end(), LessThanAbs()));
    };

    constexpr InteractionDataFrame& operator+=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i=0; i < total_size; ++i) {
            _data[i] += other._data[i];
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator-=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i=0; i < total_size; ++i) {
            _data[i] -= other._data[i];
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator*=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i=0; i < total_size; ++i) {
            _data[i] *= other._data[i];
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator/=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i=0; i < total_size; ++i) {
            _data[i] /= other._data[i];
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator*=(const double other) noexcept {
#pragma omp parallel for
        for (int i=0; i < total_size; ++i) {
            _data[i] *= other;
        }
        return *this;
    }

    constexpr InteractionDataFrame& operator/=(const double other) noexcept {
#pragma omp parallel for
        for (int i=0; i < total_size; ++i) {
            _data[i] /= other;
        }
        return *this;
    }

    constexpr std::vector<double>& get_data() noexcept {
        return _data;
    }

    constexpr const std::vector<double>& get_data() const noexcept {
        return _data;
    }

    constexpr std::vector<std::vector<std::vector<double>>> as_3D_array() const noexcept {
        std::vector<std::vector<std::vector<double>>> result(_N, std::vector<std::vector<double>>(_N, std::vector<double>(_N)));

        for (std::size_t x=0U; x < _N; ++x) {
            for (std::size_t y=0U; y < _N; ++y) {
                for (std::size_t z=0U; z < _N; ++z) {
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