#pragma once

#include "../L.hpp"
#include "../helper_functions.hpp"
#include "../momentum_iterator.hpp"
#include "MomentumGrid.hpp"
#include "MomentumSymmetry.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <omp.h>
#include <utility>
#include <vector>

namespace NickelCUT {

class InteractionDataFrame {
    using symmetry = MomentumSymmetry<N, 3>;
    static constexpr int stored_size = symmetry::independent_size();
    using storage_index_type = int;
    static_assert(stored_size <= std::numeric_limits<storage_index_type>::max());
    std::vector<double> _data;

public:
    static const std::vector<storage_index_type>& storage_indices() noexcept {
        return symmetry::storage_indices();
    }

    static std::size_t data_index(std::size_t x, std::size_t y, std::size_t z) noexcept {
        return static_cast<std::size_t>(symmetry::data_index(x, y, z));
    }

    static constexpr int independent_size() noexcept { return stored_size; }

    static const std::vector<std::size_t>& independent_indices() noexcept {
        return symmetry::independent_indices();
    }

    static const std::vector<std::size_t>& active_independent_indices() {
        static std::vector<std::size_t> indices;
        static std::uint64_t cached_generation = std::numeric_limits<std::uint64_t>::max();
        if (cached_generation != MomentumGrid::generation()) {
            indices.clear();
            for (const std::size_t index : independent_indices()) {
                const int k = static_cast<int>(index / (N * N));
                const int p = static_cast<int>((index / N) % N);
                const int q = static_cast<int>(index % N);
                if (MomentumGrid::is_fermion(k) && MomentumGrid::is_fermion(p) && MomentumGrid::is_transfer(q)) {
                    indices.push_back(index);
                }
            }
            cached_generation = MomentumGrid::generation();
        }
        return indices;
    }

    template <class Archive>
    void serialize(Archive& ar, [[maybe_unused]] const unsigned int version) {
        ar & _data;
    }

    InteractionDataFrame() : _data(stored_size, 0.0) {}
    InteractionDataFrame(double value) : _data(stored_size, value) {}

    InteractionDataFrame(const std::vector<double>& data) : _data(data) { assert(data.size() == stored_size); };
    InteractionDataFrame(std::vector<double>&& data) : _data(std::move(data)) { assert(_data.size() == stored_size); };

    InteractionDataFrame(const InteractionDataFrame&) = default;
    InteractionDataFrame& operator=(const InteractionDataFrame&) = default;

    InteractionDataFrame(InteractionDataFrame&&) noexcept = default;
    InteractionDataFrame& operator=(InteractionDataFrame&&) noexcept = default;

    inline double& operator()(std::size_t x, std::size_t y, std::size_t z) noexcept {
        assert(x < N);
        assert(y < N);
        assert(z < N);
        x = static_cast<std::size_t>(
            MomentumGrid::nearest_fermion(momentum_iterator<L>(static_cast<int>(x))).get_position());
        y = static_cast<std::size_t>(
            MomentumGrid::nearest_fermion(momentum_iterator<L>(static_cast<int>(y))).get_position());
        z = static_cast<std::size_t>(
            MomentumGrid::nearest_transfer(momentum_iterator<L>(static_cast<int>(z))).get_position());
        return _data[data_index(x, y, z)];
    }

    inline const double& operator()(std::size_t x, std::size_t y, std::size_t z) const noexcept {
        assert(x < N);
        assert(y < N);
        assert(z < N);
        x = static_cast<std::size_t>(
            MomentumGrid::nearest_fermion(momentum_iterator<L>(static_cast<int>(x))).get_position());
        y = static_cast<std::size_t>(
            MomentumGrid::nearest_fermion(momentum_iterator<L>(static_cast<int>(y))).get_position());
        z = static_cast<std::size_t>(
            MomentumGrid::nearest_transfer(momentum_iterator<L>(static_cast<int>(z))).get_position());
        return _data[data_index(x, y, z)];
    }

    inline void symmetrize() noexcept {
        for (momentum_iterator<L> K = momentum_iterator<L>::begin(); K != momentum_iterator<L>::end(); ++K) {
            for (momentum_iterator<L> P = momentum_iterator<L>::begin(); P != momentum_iterator<L>::end(); ++P) {
                for (momentum_iterator<L> Q = momentum_iterator<L>::begin(); Q != momentum_iterator<L>::end(); ++Q) {
                    operator()(K, P, Q) = 0.5 * (operator()(K, P, Q) + operator()(P, K, -Q));
                    operator()(P, K, -Q) = operator()(K, P, Q);
                }
            }
        }
    }

    inline void antisymmetrize() noexcept {
        for (momentum_iterator<L> K = momentum_iterator<L>::begin(); K != momentum_iterator<L>::end(); ++K) {
            for (momentum_iterator<L> P = momentum_iterator<L>::begin(); P != momentum_iterator<L>::end(); ++P) {
                for (momentum_iterator<L> Q = momentum_iterator<L>::begin(); Q != momentum_iterator<L>::end(); ++Q) {
                    operator()(K, P, Q) = 0.5 * (operator()(K, P, Q) - operator()(K, P, P - K - Q));
                    operator()(K, P, P - K - Q) = -operator()(K, P, Q);
                }
            }
        }
    }

    void clear_noise() noexcept {
        for (auto& val : _data) {
            if (is_zero(val))
                val = 0.;
        }
    }

    inline void fill(double value) noexcept {
        for (auto& element : _data) {
            element = value;
        }
    }

    inline double abs_squared_total() const noexcept {  ///< the square of the L2 norm
        double val{};
        for (const auto& k : MomentumGrid::fermion_points()) {
            for (const auto& p : MomentumGrid::fermion_points()) {
                for (const auto& q : MomentumGrid::transfer_points()) {
                    const double value = operator()(k.momentum, p.momentum, q.momentum);
                    val += k.weight * p.weight * q.weight * value * value;
                }
            }
        }
        return val;
    };

    inline std::size_t size() const noexcept { return _data.size(); }

    inline double abs_total() const noexcept {  ///< the L2 norm
        return std::sqrt(abs_squared_total());
    };
    inline double norm_inf() const noexcept {  ///< the L_infinity norm
        double result{};
        for (const std::size_t index : active_independent_indices()) {
            result = std::max(result, std::abs(_data[data_index(
                index / (N * N), (index / N) % N, index % N)]));
        }
        return result;
    }

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
        for (int i = 0; i < stored_size; ++i) {
            _data[i] += other._data[i];
        }
        return *this;
    }

    inline InteractionDataFrame& operator-=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i = 0; i < stored_size; ++i) {
            _data[i] -= other._data[i];
        }
        return *this;
    }

    inline InteractionDataFrame& operator*=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i = 0; i < stored_size; ++i) {
            _data[i] *= other._data[i];
        }
        return *this;
    }

    inline InteractionDataFrame& operator/=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i = 0; i < stored_size; ++i) {
            _data[i] /= other._data[i];
        }
        return *this;
    }

    inline InteractionDataFrame& operator*=(const double other) noexcept {
#pragma omp parallel for
        for (int i = 0; i < stored_size; ++i) {
            _data[i] *= other;
        }
        return *this;
    }

    inline InteractionDataFrame& operator/=(const double other) noexcept {
#pragma omp parallel for
        for (int i = 0; i < stored_size; ++i) {
            _data[i] /= other;
        }
        return *this;
    }

    inline InteractionDataFrame& operator+=(const double other) noexcept {
#pragma omp parallel for
        for (int i = 0; i < stored_size; ++i) {
            _data[i] += other;
        }
        return *this;
    }

    inline std::vector<double>& get_data() noexcept { return _data; }

    inline const std::vector<double>& get_data() const noexcept { return _data; }

    inline std::vector<std::vector<std::vector<double>>> as_3D_array() const noexcept {
        std::vector<std::vector<std::vector<double>>> result(
            N, std::vector<std::vector<double>>(N, std::vector<double>(N)));

        for (std::size_t x = 0U; x < N; ++x) {
            for (std::size_t y = 0U; y < N; ++y) {
                for (std::size_t z = 0U; z < N; ++z) {
                    result[x][y][z] = (*this)(x, y, z);
                }
            }
        }
        return result;
    }
};

template <int N>
inline InteractionDataFrame operator+(InteractionDataFrame lhs, const InteractionDataFrame& rhs) {
    return (lhs += rhs);
}
template <int N>
inline InteractionDataFrame operator-(InteractionDataFrame lhs, const InteractionDataFrame& rhs) {
    return (lhs -= rhs);
}
template <int N>
inline InteractionDataFrame operator*(InteractionDataFrame lhs, const InteractionDataFrame& rhs) {
    return (lhs *= rhs);
}
template <int N>
inline InteractionDataFrame operator/(InteractionDataFrame lhs, const InteractionDataFrame& rhs) {
    return (lhs /= rhs);
}

template <int N>
inline InteractionDataFrame operator*(InteractionDataFrame lhs, const double rhs) {
    return (lhs *= rhs);
}
template <int N>
inline InteractionDataFrame operator*(const double lhs, InteractionDataFrame rhs) {
    return (rhs *= lhs);
}
template <int N>
inline InteractionDataFrame operator/(InteractionDataFrame lhs, const double rhs) {
    return (lhs /= rhs);
}
template <int N>
inline InteractionDataFrame operator+(InteractionDataFrame lhs, const double rhs) {
    return (lhs += rhs);
}
template <int N>
inline InteractionDataFrame operator+(const double lhs, InteractionDataFrame rhs) {
    return (rhs += lhs);
}

}  // namespace NickelCUT