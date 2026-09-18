#pragma once

#include "../L.hpp"
#include "InteractionDataFrame.hpp"
#include "momentum_iterator.hpp"

#include <nlohmann/json.hpp>

#include <array>

namespace NickelCUT::flow {

struct Model;

struct FlowContainer {
    using coeff_array = std::array<double, N>;
    using mom_it = NickelCUT::flow::momentum_iterator<L>;

    InteractionDataFrame interactions_same_spin;
    InteractionDataFrame interactions_differing_spin;
    coeff_array dispersion;
    coeff_array epsilon_tilde;

    template<class Archive>
    void serialize(Archive& ar, [[maybe_unused]] const unsigned int version) {
        ar & interactions_same_spin;
        ar & interactions_differing_spin;
        ar & dispersion;
        ar & epsilon_tilde;
    }

    FlowContainer() noexcept;
    // Note that the way we write down our model (sum_(sigma) U c_(sigma)^dagger c_(-sigma)^dagger c_(-sigma) c_(sigma))
    // Introduces an additional factor of 2 compared to the standard Hubbard model, in which there is no spin summation
    // in the interaction term. This additional factor of 2 is explicitly removed when invoking this constructor.
    // That is, if model.U_0 = 1, this function will fill the InteractionDataFrame with 0.5
    // Thereby, we can use the same notation as in the standard Hubbard model without introducing unnecessary confusion.
    FlowContainer(const Model& model) noexcept;

    FlowContainer(const FlowContainer&) = default;
    FlowContainer& operator=(const FlowContainer&) = default;

    FlowContainer(FlowContainer&&) noexcept = default;
    FlowContainer& operator=(FlowContainer&&) noexcept = default;

    void fill_epsilon_tilde() noexcept;

    void fill(double value) noexcept;
    void reset() noexcept;
    bool contains_nan_or_inf() const noexcept;
    double residual_offdiagonality() const noexcept;

    bool is_inversion_symmetric() const noexcept;
    bool is_hermitian() const noexcept;
    bool is_particle_hole_invariant() const noexcept;

    double abs_total() const noexcept; ///< the L2 norm
    double norm_inf() const noexcept; ///< the L_infinity norm

    // Required for boost odeint
    FlowContainer& operator+=(const FlowContainer& other) noexcept;
    FlowContainer& operator-=(const FlowContainer& other) noexcept;
    FlowContainer& operator*=(const FlowContainer& other) noexcept;
    FlowContainer& operator/=(const FlowContainer& other) noexcept;

    FlowContainer& operator*=(const double other) noexcept;
    FlowContainer& operator/=(const double other) noexcept;
    FlowContainer& operator+=(const double other) noexcept;

private:
    void print_memory_allocation() const;
};

inline FlowContainer operator+(FlowContainer lhs, const FlowContainer& rhs) noexcept { return (lhs += rhs); }
inline FlowContainer operator-(FlowContainer lhs, const FlowContainer& rhs) noexcept { return (lhs -= rhs); }
inline FlowContainer operator*(FlowContainer lhs, const FlowContainer& rhs) noexcept { return (lhs *= rhs); }
inline FlowContainer operator/(FlowContainer lhs, const FlowContainer& rhs) noexcept { return (lhs /= rhs); }

inline FlowContainer operator*(FlowContainer lhs, const double rhs) noexcept { return (lhs *= rhs); }
inline FlowContainer operator*(const double lhs, FlowContainer rhs) noexcept { return (rhs *= lhs); }
inline FlowContainer operator/(FlowContainer lhs, const double rhs) noexcept { return (lhs /= rhs); }
inline FlowContainer operator+(FlowContainer lhs, const double rhs) noexcept { return (lhs += rhs); }
inline FlowContainer operator+(const double lhs, FlowContainer rhs) noexcept { return (rhs += lhs); }

FlowContainer abs(FlowContainer input) noexcept; ///< elementwise L2 norm

void to_json(nlohmann::json& j, const FlowContainer& container) noexcept;

} // namespace NickelCUT::flow