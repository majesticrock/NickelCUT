#pragma once

#include "L.hpp"
#include "FlowContainerIterator.hpp"
#include "InteractionDataFrame.hpp"
#include "momentum_iterator.hpp"
#include "Model.hpp"

#include <nlohmann/json.hpp>

#include <array>
#include <iostream>

namespace NickelCUT::flow {

struct FlowContainer {
    using coeff_array = std::array<double, N>;
    using mom_it = NickelCUT::flow::momentum_iterator<L>;

    InteractionDataFrame<N> interactions_same_spin;
    InteractionDataFrame<N> interactions_differing_spin;
    coeff_array dispersion;
    coeff_array epsilon_tilde;
    coeff_array occupation_numbers;

    FlowContainer();
    FlowContainer(const Model& model);

    void fill_epsilon_tilde();

    void fill(double value) noexcept;
    void reset() noexcept;
    bool contains_nan_or_inf() const noexcept;

    // Required for boost odeint
    double abs() const; ///< uses the L1 norm
    double norm_inf() const; ///< uses the L_infinity norm

    FlowContainer& operator+=(const FlowContainer& other);
    FlowContainer& operator-=(const FlowContainer& other);
    FlowContainer& operator*=(const FlowContainer& other);
    FlowContainer& operator/=(const FlowContainer& other);

    FlowContainer& operator*=(const double other);
    FlowContainer& operator/=(const double other);

    using iterator = FlowContainerIterator<FlowContainer>;
    using const_iterator = FlowContainerIterator<const FlowContainer>;

    constexpr iterator begin() noexcept { return iterator(this, 0); }
    constexpr iterator end() noexcept { return iterator(this, total_value_count); }

    constexpr const_iterator begin() const noexcept { return const_iterator(this, 0); }
    constexpr const_iterator end() const noexcept { return const_iterator(this, total_value_count); }
    constexpr const_iterator cbegin() const noexcept { return begin(); }
    constexpr const_iterator cend() const noexcept { return end(); }

    constexpr double& value_at(std::size_t idx) noexcept {
        if (idx < interaction_size) {
            return interactions_same_spin.get_data()[idx];
        }
        idx -= interaction_size;
        if (idx < interaction_size) {
            return interactions_differing_spin.get_data()[idx];
        }
        idx -= interaction_size;
        return dispersion[idx];
    }

    constexpr const double& value_at(std::size_t idx) const noexcept {
        if (idx < interaction_size) {
            return interactions_same_spin.get_data()[idx];
        }
        idx -= interaction_size;
        if (idx < interaction_size) {
            return interactions_differing_spin.get_data()[idx];
        }
        idx -= interaction_size;
        return dispersion[idx];
    }

    static constexpr std::size_t interaction_size = N * N * N;
    static constexpr std::size_t total_value_count = 2 * interaction_size + N;
};

inline FlowContainer operator+(FlowContainer lhs, const FlowContainer& rhs) { return (lhs += rhs); }
inline FlowContainer operator-(FlowContainer lhs, const FlowContainer& rhs) { return (lhs -= rhs); }
inline FlowContainer operator*(FlowContainer lhs, const FlowContainer& rhs) { return (lhs *= rhs); }
inline FlowContainer operator/(FlowContainer lhs, const FlowContainer& rhs) { return (lhs /= rhs); }

inline FlowContainer operator*(FlowContainer lhs, const double rhs) { return (lhs *= rhs); }
inline FlowContainer operator*(const double lhs, FlowContainer rhs) { return (rhs *= lhs); }
inline FlowContainer operator/(FlowContainer lhs, const double rhs) { return (lhs /= rhs); }

void to_json(nlohmann::json& j, const FlowContainer& container);

} // namespace NickelCUT::flow