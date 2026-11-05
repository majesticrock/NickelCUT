#pragma once

#include "../helper_functions.hpp"
#include "../momentum_iterator.hpp"
#include "MomentumGrid.hpp"
#include "MomentumSymmetry.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <omp.h>
#include <utility>
#include <vector>

namespace NickelCUT {

class SingleParticleDataFrame {
    using symmetry = MomentumSymmetry<N, 1>;
    static constexpr int stored_size = symmetry::independent_size();
    std::vector<double> _data;

public:
    static const std::vector<int>& storage_indices() noexcept { return symmetry::storage_indices(); }

    static constexpr int independent_size() noexcept { return stored_size; }

    static const std::vector<std::size_t>& independent_indices() noexcept { return symmetry::independent_indices(); }

    static const std::vector<std::size_t>& active_independent_indices() {
        static std::vector<std::size_t> indices;
        static std::uint64_t cached_generation = std::numeric_limits<std::uint64_t>::max();
        if (cached_generation != MomentumGrid::generation()) {
            indices.clear();
            for (const std::size_t index : independent_indices()) {
                if (MomentumGrid::is_fermion(static_cast<int>(index))) {
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

    SingleParticleDataFrame() : _data(stored_size, 0.0) {}
    SingleParticleDataFrame(double value) : _data(stored_size, value) {}
    SingleParticleDataFrame(const std::vector<double>& data) : _data(data) { assert(data.size() == stored_size); }
    SingleParticleDataFrame(std::vector<double>&& data) : _data(std::move(data)) {
        assert(_data.size() == stored_size);
    }

    SingleParticleDataFrame(const SingleParticleDataFrame&) = default;
    SingleParticleDataFrame& operator=(const SingleParticleDataFrame&) = default;
    SingleParticleDataFrame(SingleParticleDataFrame&&) noexcept = default;
    SingleParticleDataFrame& operator=(SingleParticleDataFrame&&) noexcept = default;

    double& operator()(std::size_t momentum) noexcept {
        assert(momentum < N);
        momentum = static_cast<std::size_t>(
            MomentumGrid::nearest_fermion(momentum_iterator<L>(static_cast<int>(momentum))).get_position());
        return _data[symmetry::data_index(momentum)];
    }

    const double& operator()(std::size_t momentum) const noexcept {
        assert(momentum < N);
        momentum = static_cast<std::size_t>(
            MomentumGrid::nearest_fermion(momentum_iterator<L>(static_cast<int>(momentum))).get_position());
        return _data[symmetry::data_index(momentum)];
    }

    double& operator[](std::size_t momentum) noexcept { return (*this)(momentum); }
    const double& operator[](std::size_t momentum) const noexcept { return (*this)(momentum); }

    void clear_noise() noexcept {
        for (auto& value : _data) {
            if (is_zero(value)) {
                value = 0.0;
            }
        }
    }

    void fill(double value) noexcept { std::fill(_data.begin(), _data.end(), value); }

    double abs_squared_total() const noexcept {
        double result = 0.0;
        for (const auto& point : MomentumGrid::fermion_points()) {
            const double value = (*this)(point.momentum);
            result += point.weight * value * value;
        }
        return result;
    }

    std::size_t size() const noexcept { return _data.size(); }
    double abs_total() const noexcept { return std::sqrt(abs_squared_total()); }
    double norm_inf() const noexcept {
        double result{};
        for (const auto& point : MomentumGrid::fermion_points()) {
            result = std::max(result, std::abs((*this)(point.momentum)));
        }
        return result;
    }

    SingleParticleDataFrame abs() const noexcept {
        SingleParticleDataFrame result(*this);
        for (auto& value : result._data) {
            value = std::abs(value);
        }
        return result;
    }

    SingleParticleDataFrame& abs_in_place() noexcept {
        for (auto& value : _data) {
            value = std::abs(value);
        }
        return *this;
    }

    SingleParticleDataFrame& operator+=(const SingleParticleDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i = 0; i < stored_size; ++i) {
            _data[i] += other._data[i];
        }
        return *this;
    }

    SingleParticleDataFrame& operator-=(const SingleParticleDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i = 0; i < stored_size; ++i) {
            _data[i] -= other._data[i];
        }
        return *this;
    }

    SingleParticleDataFrame& operator*=(const SingleParticleDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i = 0; i < stored_size; ++i) {
            _data[i] *= other._data[i];
        }
        return *this;
    }

    SingleParticleDataFrame& operator/=(const SingleParticleDataFrame& other) noexcept {
#pragma omp parallel for
        for (int i = 0; i < stored_size; ++i) {
            _data[i] /= other._data[i];
        }
        return *this;
    }

    SingleParticleDataFrame& operator*=(double other) noexcept {
#pragma omp parallel for
        for (int i = 0; i < stored_size; ++i) {
            _data[i] *= other;
        }
        return *this;
    }

    SingleParticleDataFrame& operator/=(double other) noexcept {
#pragma omp parallel for
        for (int i = 0; i < stored_size; ++i) {
            _data[i] /= other;
        }
        return *this;
    }

    SingleParticleDataFrame& operator+=(double other) noexcept {
#pragma omp parallel for
        for (int i = 0; i < stored_size; ++i) {
            _data[i] += other;
        }
        return *this;
    }

    std::vector<double>& get_data() noexcept { return _data; }
    const std::vector<double>& get_data() const noexcept { return _data; }

    std::array<double, N> as_array() const {
        std::array<double, N> result{};
        for (std::size_t momentum = 0; momentum < N; ++momentum) {
            result[momentum] = (*this)(momentum);
        }
        return result;
    }
};

inline SingleParticleDataFrame operator+(SingleParticleDataFrame lhs, const SingleParticleDataFrame& rhs) {
    return lhs += rhs;
}
inline SingleParticleDataFrame operator-(SingleParticleDataFrame lhs, const SingleParticleDataFrame& rhs) {
    return lhs -= rhs;
}
inline SingleParticleDataFrame operator*(SingleParticleDataFrame lhs, const SingleParticleDataFrame& rhs) {
    return lhs *= rhs;
}
inline SingleParticleDataFrame operator/(SingleParticleDataFrame lhs, const SingleParticleDataFrame& rhs) {
    return lhs /= rhs;
}
inline SingleParticleDataFrame operator*(SingleParticleDataFrame lhs, double rhs) { return lhs *= rhs; }
inline SingleParticleDataFrame operator*(double lhs, SingleParticleDataFrame rhs) { return rhs *= lhs; }
inline SingleParticleDataFrame operator/(SingleParticleDataFrame lhs, double rhs) { return lhs /= rhs; }
inline SingleParticleDataFrame operator+(SingleParticleDataFrame lhs, double rhs) { return lhs += rhs; }
inline SingleParticleDataFrame operator+(double lhs, SingleParticleDataFrame rhs) { return rhs += lhs; }

}  // namespace NickelCUT
