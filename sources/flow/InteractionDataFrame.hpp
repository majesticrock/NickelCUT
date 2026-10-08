#pragma once

#include "../helper_functions.hpp"
#include "momentum_iterator.hpp"
#include "../L.hpp"

#include <omp.h>

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <algorithm>
#include <vector>
#include <iostream>

namespace NickelCUT::flow {

class InteractionDataFrame {
    static constexpr int total_size = N * N * N;
    static constexpr int fixed_axis_reflection_size = (2 * L) * (2 * L) * (2 * L);
    static constexpr int fixed_both_reflections_size = 4 * 4 * 4;
    static constexpr int stored_size = (total_size + 2 * fixed_axis_reflection_size + fixed_both_reflections_size) / 4;
    using storage_index_type = int;
    static_assert(stored_size <= std::numeric_limits<storage_index_type>::max());
    std::vector<double> _data;

    static std::size_t reflect_kx(std::size_t momentum) noexcept {
        const std::size_t kx = momentum % L;
        const std::size_t ky = momentum / L;
        return (static_cast<std::size_t>(L) - kx) % L + L * ky;
    }

    static std::size_t reflect_ky(std::size_t momentum) noexcept {
        const std::size_t kx = momentum % L;
        const std::size_t ky = momentum / L;
        return kx + L * ((static_cast<std::size_t>(L) - ky) % L);
    }

    static std::size_t canonical_index(std::size_t index) noexcept {
        const std::size_t x = index / (N * N);
        const std::size_t y = (index / N) % N;
        const std::size_t z = index % N;
        const auto flatten = [](std::size_t a, std::size_t b, std::size_t c) {
            return c + N * (b + N * a);
        };

        return std::min({
            index,
            flatten(reflect_kx(x), reflect_kx(y), reflect_kx(z)),
            flatten(reflect_ky(x), reflect_ky(y), reflect_ky(z)),
            flatten(reflect_kx(reflect_ky(x)), reflect_kx(reflect_ky(y)), reflect_kx(reflect_ky(z)))
        });
    }

    static const std::array<storage_index_type, total_size>& storage_indices() noexcept {
        static const std::array<storage_index_type, total_size> indices = [] {
            std::array<storage_index_type, total_size> result{};
            std::array<storage_index_type, total_size> canonical_to_storage{};
            canonical_to_storage.fill(std::numeric_limits<storage_index_type>::max());

            std::size_t next_storage_index = 0;
            for (std::size_t index = 0; index < total_size; ++index) {
                const std::size_t canonical = canonical_index(index);
                if (canonical == index) {
                    canonical_to_storage[canonical] = static_cast<storage_index_type>(next_storage_index++);
                }
            }
            assert(next_storage_index == stored_size);

            for (std::size_t index = 0; index < total_size; ++index) {
                result[index] = canonical_to_storage[canonical_index(index)];
            }
            return result;
        }();
        return indices;
    }

    static std::size_t data_index(std::size_t x, std::size_t y, std::size_t z) noexcept {
        return storage_indices()[z + N * (y + N * x)];
    }

public:
    static constexpr std::size_t independent_size() noexcept {
        return stored_size;
    }

    static const std::array<std::size_t, stored_size>& independent_indices() noexcept {
        static const std::array<std::size_t, stored_size> indices = [] {
            std::array<std::size_t, stored_size> result{};
            std::size_t next_index = 0;
            for (std::size_t index = 0; index < total_size; ++index) {
                if (canonical_index(index) == index) {
                    result[next_index++] = index;
                }
            }
            assert(next_index == stored_size);
            return result;
        }();
        return indices;
    }

    template<class Archive>
    void serialize(Archive& ar, [[maybe_unused]] const unsigned int version) {
        ar & _data;
    }

    InteractionDataFrame()
        : _data(stored_size, 0.0)
    {}
    InteractionDataFrame(double value)
        : _data(stored_size, value)
    {}

    InteractionDataFrame(const std::vector<double>& data)
        : _data(data) { assert(data.size() == stored_size); };
    InteractionDataFrame(std::vector<double>&& data)
        : _data(std::move(data)) { assert(_data.size() == stored_size); };

    InteractionDataFrame(const InteractionDataFrame&) = default;
    InteractionDataFrame& operator=(const InteractionDataFrame&) = default;
    
    InteractionDataFrame(InteractionDataFrame&&) noexcept = default;
    InteractionDataFrame& operator=(InteractionDataFrame&&) noexcept = default;

    inline double& operator()(std::size_t x,
                                 std::size_t y,
                                 std::size_t z) noexcept {
        assert(x < N);
        assert(y < N);
        assert(z < N);
        return _data[data_index(x, y, z)];
    }

    inline const double& operator()(std::size_t x,
                                       std::size_t y,
                                       std::size_t z) const noexcept {
        assert(x < N);
        assert(y < N);
        assert(z < N);
        return _data[data_index(x, y, z)];
    }

    inline void symmetrize() noexcept {
        for (momentum_iterator<L> K = momentum_iterator<L>::begin(); K != momentum_iterator<L>::end(); ++K) {
        for (momentum_iterator<L> P = momentum_iterator<L>::begin(); P != momentum_iterator<L>::end(); ++P) {
        for (momentum_iterator<L> Q = momentum_iterator<L>::begin(); Q != momentum_iterator<L>::end(); ++Q) {
            operator()(K, P, Q) = 0.5 * (
                operator()(K, P, Q) + operator()(P, K, -Q)
            );
            operator()(P, K, -Q) = operator()(K, P, Q);
        }}}
    }

    inline void antisymmetrize() noexcept {
        for (momentum_iterator<L> K = momentum_iterator<L>::begin(); K != momentum_iterator<L>::end(); ++K) {
        for (momentum_iterator<L> P = momentum_iterator<L>::begin(); P != momentum_iterator<L>::end(); ++P) {
        for (momentum_iterator<L> Q = momentum_iterator<L>::begin(); Q != momentum_iterator<L>::end(); ++Q) {
            operator()(K, P, Q) = 0.5 * (
                operator()(K, P, Q) - operator()(K, P, P-K-Q)
            );
            operator()(K, P, P-K-Q) = -operator()(K, P, Q);
        }}}
    }

    void clear_noise() noexcept {
        for (auto& val : _data) {
            if (is_zero(val)) val = 0.;
        }
    }

    inline void fill(double value) noexcept {
        for (auto& element : _data) {
            element = value;
        }
    }

    inline double abs_squared_total() const noexcept {  ///< the square of the L2 norm
        double val{};
        const auto& indices = storage_indices();
        for (const std::size_t index : indices) {
            val += _data[index] * _data[index];
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
        for (int i=0; i < stored_size; ++i) {
            _data[i] += other._data[i];
        }
        return *this;
    }

    inline InteractionDataFrame& operator-=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i=0; i < stored_size; ++i) {
            _data[i] -= other._data[i];
        }
        return *this;
    }

    inline InteractionDataFrame& operator*=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i=0; i < stored_size; ++i) {
            _data[i] *= other._data[i];
        }
        return *this;
    }

    inline InteractionDataFrame& operator/=(const InteractionDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i=0; i < stored_size; ++i) {
            _data[i] /= other._data[i];
        }
        return *this;
    }

    inline InteractionDataFrame& operator*=(const double other) noexcept {
#pragma omp parallel for
        for (int i=0; i < stored_size; ++i) {
            _data[i] *= other;
        }
        return *this;
    }

    inline InteractionDataFrame& operator/=(const double other) noexcept {
#pragma omp parallel for
        for (int i=0; i < stored_size; ++i) {
            _data[i] /= other;
        }
        return *this;
    }

    inline InteractionDataFrame& operator+=(const double other) noexcept {
#pragma omp parallel for
        for (int i=0; i < stored_size; ++i) {
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