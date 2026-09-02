#pragma once

#include "../helper_functions.hpp"
#include "momentum_iterator.hpp"
#include "../L.hpp"

#include <omp.h>

#include <cassert>
#include <cstddef>
#include <utility>
#include <algorithm>
#include <vector>
#include <iostream>

namespace NickelCUT::flow {

class InteractionDataFrame {
    static constexpr int total_size = N * N * N;
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

    InteractionDataFrame(const std::vector<double>& data)
        : _data(data) { assert(data.size() == total_size); };
    InteractionDataFrame(std::vector<double>&& data)
        : _data(std::move(data)) { assert(data.size() == total_size); };

    inline double& operator()(std::size_t x,
                                 std::size_t y,
                                 std::size_t z) noexcept {
        assert(x < N);
        assert(y < N);
        assert(z < N);
        return _data[z + N * (y + N * x)];
    }

    inline const double& operator()(std::size_t x,
                                       std::size_t y,
                                       std::size_t z) const noexcept {
        assert(x < N);
        assert(y < N);
        assert(z < N);
        return _data[z + N * (y + N * x)];
    }

    inline void symmetrize() noexcept {
        for (momentum_iterator<L> K = momentum_iterator<L>::begin(); K != momentum_iterator<L>::end(); ++K) {
        for (momentum_iterator<L> P = momentum_iterator<L>::begin(); P != momentum_iterator<L>::end(); ++P) {
        for (momentum_iterator<L> Q = momentum_iterator<L>::begin(); Q != momentum_iterator<L>::end(); ++Q) {
            operator()(K.get_position(), P.get_position(), Q.get_position()) = 0.5 * (
                operator()(K.get_position(), P.get_position(), Q.get_position()) 
                + operator()(P.get_position(), K.get_position(), (-Q).get_position())
            );
            operator()(P.get_position(), K.get_position(), (-Q).get_position()) = operator()(K.get_position(), P.get_position(), Q.get_position());
        }}}
    }

    inline void fill(double value) noexcept {
        for (auto& element : _data) {
            element = value;
        }
    }

    inline double abs_squared_total() const noexcept {  ///< the square of the L2 norm
        double val{};
        for(const auto& element : _data) {
            val += element*element;
        }
        return val;
    };

    inline std::size_t size() const noexcept {
        return _data.size();
    }

    inline double abs_total() const noexcept {  ///< the L2 norm
        return std::sqrt(abs_squared_total());
    };
    inline double norm_inf() const noexcept {  ///< the L_infinity norm
        return std::abs(*std::max_element(_data.begin(), _data.end(), LessThanAbs()));
    };

    // Required for boost odeint
    inline InteractionDataFrame abs() const noexcept {
        InteractionDataFrame ret(*this);
        for (auto& val : ret._data) {
            val = std::abs(val);
        }
        return ret;
    }
    inline InteractionDataFrame& abs_in_place() noexcept {
        for (auto& val : _data) {
            val = std::abs(val);
        }
        return *this;
    }

    inline InteractionDataFrame& operator+=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i=0; i < total_size; ++i) {
            _data[i] += other._data[i];
        }
        return *this;
    }

    inline InteractionDataFrame& operator-=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i=0; i < total_size; ++i) {
            _data[i] -= other._data[i];
        }
        return *this;
    }

    inline InteractionDataFrame& operator*=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i=0; i < total_size; ++i) {
            _data[i] *= other._data[i];
        }
        return *this;
    }

    inline InteractionDataFrame& operator/=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i=0; i < total_size; ++i) {
            _data[i] /= other._data[i];
        }
        return *this;
    }

    inline InteractionDataFrame& operator*=(const double other) noexcept {
#pragma omp parallel for
        for (int i=0; i < total_size; ++i) {
            _data[i] *= other;
        }
        return *this;
    }

    inline InteractionDataFrame& operator/=(const double other) noexcept {
#pragma omp parallel for
        for (int i=0; i < total_size; ++i) {
            _data[i] /= other;
        }
        return *this;
    }

    inline InteractionDataFrame& operator+=(const double other) noexcept {
#pragma omp parallel for
        for (int i=0; i < total_size; ++i) {
            _data[i] += other;
        }
        return *this;
    }

    inline std::vector<double>& get_data() noexcept {
        return _data;
    }

    inline const std::vector<double>& get_data() const noexcept {
        return _data;
    }

    inline std::vector<std::vector<std::vector<double>>> as_3D_array() const noexcept {
        std::vector<std::vector<std::vector<double>>> result(N, std::vector<std::vector<double>>(N, std::vector<double>(N)));

        for (std::size_t x=0U; x < N; ++x) {
            for (std::size_t y=0U; y < N; ++y) {
                for (std::size_t z=0U; z < N; ++z) {
                    result[x][y][z] = (*this)(x, y, z);
                }
            }
        }
        return result;
    }
};

template <int N>
inline InteractionDataFrame operator+(InteractionDataFrame lhs, const InteractionDataFrame& rhs) { return (lhs += rhs); }
template <int N>
inline InteractionDataFrame operator-(InteractionDataFrame lhs, const InteractionDataFrame& rhs) { return (lhs -= rhs); }
template <int N>
inline InteractionDataFrame operator*(InteractionDataFrame lhs, const InteractionDataFrame& rhs) { return (lhs *= rhs); }
template <int N>
inline InteractionDataFrame operator/(InteractionDataFrame lhs, const InteractionDataFrame& rhs) { return (lhs /= rhs); }

template <int N>
inline InteractionDataFrame operator*(InteractionDataFrame lhs, const double rhs) { return (lhs *= rhs); }
template <int N>
inline InteractionDataFrame operator*(const double lhs, InteractionDataFrame rhs) { return (rhs *= lhs); }
template <int N>
inline InteractionDataFrame operator/(InteractionDataFrame lhs, const double rhs) { return (lhs /= rhs); }
template <int N>
inline InteractionDataFrame operator+(InteractionDataFrame lhs, const double rhs) { return (lhs += rhs); }
template <int N>
inline InteractionDataFrame operator+(const double lhs, InteractionDataFrame rhs) { return (rhs += lhs); }

} // namespace NickelCut::flow