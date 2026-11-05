#include "FlowContainer.hpp"

#include "../helper_functions.hpp"
#include "Model.hpp"
#include "MomentumGrid.hpp"
#include "../momentum_iterator.hpp"
#include "occupation_numbers.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <iostream>
#include <limits>

namespace NickelCUT {

FlowContainer::FlowContainer() noexcept : interactions_same_spin(), interactions_differing_spin() {
    dispersion.fill(double{});
}

FlowContainer::FlowContainer(const Model& model) noexcept
    : interactions_same_spin(),
      // Factor 1/2 to account for the extra spin summation compared to the standard Hubbard model
      interactions_differing_spin(0.5 * model.U_0 / N) {
    for (const auto& point : MomentumGrid::fermion_points()) {
        const auto K = point.momentum;
        dispersion[K] = model.epsilon_0(K.get_kx(), K.get_ky());
    }
}

void FlowContainer::print_memory_allocation() const {
    const std::size_t allocated_bytes = sizeof(FlowContainer) +
                                        interactions_same_spin.get_data().capacity() * sizeof(double) +
                                        interactions_differing_spin.get_data().capacity() * sizeof(double) +
                                        dispersion.get_data().capacity() * sizeof(double);

    std::cout << "FlowContainer constructed, allocated memory: " << allocated_bytes / (1024) << " kb" << std::endl;
}

void FlowContainer::fill(double value) noexcept {
    interactions_same_spin.fill(value);
    interactions_differing_spin.fill(value);
    dispersion.fill(value);
}

void FlowContainer::reset() noexcept {
    fill(double{});
}

bool FlowContainer::contains_nan_or_inf() const noexcept {
    for (const auto& val : interactions_same_spin.get_data()) {
        if (std::isnan(val) || std::isinf(val))
            return true;
    }
    for (const auto& val : interactions_differing_spin.get_data()) {
        if (std::isnan(val) || std::isinf(val))
            return true;
    }
    for (const auto& val : dispersion.get_data()) {
        if (std::isnan(val) || std::isinf(val))
            return true;
    }
    return false;
}

double FlowContainer::residual_offdiagonality() const noexcept {
    // A factor of 1/N is contained in the values of the interactions
    // so -> abs_squared is 1/N^2, there a 3 sums of N elements -> scales as N
    // sqrt(N) = L -> dividing by L makes ROD independent of the system size
    double val{};
    for (const auto& K_point : MomentumGrid::fermion_points()) {
        const auto K = K_point.momentum;
        for (const auto& P_point : MomentumGrid::fermion_points()) {
            const auto P = P_point.momentum;
            for (const auto& Q_point : MomentumGrid::transfer_points()) {
                const auto Q = Q_point.momentum;
                const double alpha = alpha_factor(K, P, Q);
                val += K_point.weight * P_point.weight * Q_point.weight * alpha * alpha *
                       (interactions_differing_spin(K, P, Q) * interactions_differing_spin(K, P, Q) +
                        interactions_same_spin(K, P, Q) * interactions_same_spin(K, P, Q));
            }
        }
    }
    return std::sqrt(val) / L;
}

double FlowContainer::max_generator_element() const noexcept {
    double max = 0.;
    for (const std::size_t tuple_index : InteractionDataFrame::active_independent_indices()) {
        momentum_iterator<L> K(static_cast<int>(tuple_index / (N * N)));
        momentum_iterator<L> P(static_cast<int>((tuple_index / N) % N));
        momentum_iterator<L> Q(static_cast<int>(tuple_index % N));
        double val = alpha_factor(K, P, Q) * std::max(
            std::abs(interactions_differing_spin(K,P,Q)),
            std::abs(interactions_same_spin(K,P,Q))
        );
        if (val > max) max = val;
    }
    return max;
}

bool FlowContainer::is_inversion_symmetric() const noexcept {
    for (const auto& point : MomentumGrid::fermion_points()) {
        const auto K = point.momentum;
        if (!float_equal(dispersion[K], dispersion[-K])) {
            std::cerr << "Dispersion is not inversion symmetric!" << std::endl;
            return false;
        }
    }

    for (const auto& K_point : MomentumGrid::fermion_points()) {
        const auto K = K_point.momentum;
        for (const auto& P_point : MomentumGrid::fermion_points()) {
            const auto P = P_point.momentum;
            for (const auto& Q_point : MomentumGrid::transfer_points()) {
                const auto Q = Q_point.momentum;
                if (!float_equal(interactions_differing_spin(K, P, Q), interactions_differing_spin(-K, -P, -Q))) {
                    std::cerr << "Interaction antiparallel is not inversion symmetric!" << std::endl;
                    return false;
                }
            }
        }
    }

    for (const auto& K_point : MomentumGrid::fermion_points()) {
        const auto K = K_point.momentum;
        for (const auto& P_point : MomentumGrid::fermion_points()) {
            const auto P = P_point.momentum;
            for (const auto& Q_point : MomentumGrid::transfer_points()) {
                const auto Q = Q_point.momentum;
                if (!float_equal(interactions_same_spin(K, P, Q), interactions_same_spin(-K, -P, -Q))) {
                    std::cerr << "Interaction parallel is not inversion symmetric!" << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}

bool FlowContainer::is_hermitian() const noexcept {
    for (const auto& K_point : MomentumGrid::fermion_points()) {
        const auto K = K_point.momentum;
        for (const auto& P_point : MomentumGrid::fermion_points()) {
            const auto P = P_point.momentum;
            for (const auto& Q_point : MomentumGrid::transfer_points()) {
                const auto Q = Q_point.momentum;
                if (!float_equal(interactions_differing_spin(K, P, Q), interactions_differing_spin(K + Q, P - Q, -Q))) {
                    std::cerr << "Interaction antiparallel is not Hermitian!" << std::endl;
                    return false;
                }
            }
        }
    }

    for (const auto& K_point : MomentumGrid::fermion_points()) {
        const auto K = K_point.momentum;
        for (const auto& P_point : MomentumGrid::fermion_points()) {
            const auto P = P_point.momentum;
            for (const auto& Q_point : MomentumGrid::transfer_points()) {
                const auto Q = Q_point.momentum;
                if (!float_equal(interactions_same_spin(K, P, Q), interactions_same_spin(K + Q, P - Q, -Q))) {
                    std::cerr << "Interaction parallel is not Hermitian!" << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}

bool FlowContainer::is_particle_hole_invariant() const noexcept {
    for (const auto& point : MomentumGrid::fermion_points()) {
        const auto K = point.momentum;
        if (!is_zero(dispersion[K] + dispersion[PI<L> - K])) {
            std::cerr << "Dispersion is not particle-hole invariant!       " << K << ":\t" << dispersion[K] << "    "
                      << dispersion[PI<L> - K] << std::endl;
            return false;
        }
    }

    for (const auto& K_point : MomentumGrid::fermion_points()) {
        const auto K = K_point.momentum;
        for (const auto& P_point : MomentumGrid::fermion_points()) {
            const auto P = P_point.momentum;
            for (const auto& Q_point : MomentumGrid::transfer_points()) {
                const auto Q = Q_point.momentum;
                if (!float_equal(interactions_differing_spin(K, P, Q),
                                 interactions_differing_spin(PI<L> - K - Q, PI<L> - P + Q, Q))) {
                    std::cerr << "Interaction antiparallel is not particle-hole invariant!   " << K << " " << P << " "
                              << Q << ":    " << interactions_differing_spin(K, P, Q) << "   "
                              << interactions_differing_spin(PI<L> - K - Q, PI<L> - P + Q, Q) << std::endl;

                    return false;
                }
            }
        }
    }

    for (const auto& K_point : MomentumGrid::fermion_points()) {
        const auto K = K_point.momentum;
        for (const auto& P_point : MomentumGrid::fermion_points()) {
            const auto P = P_point.momentum;
            for (const auto& Q_point : MomentumGrid::transfer_points()) {
                const auto Q = Q_point.momentum;
                if (!float_equal(interactions_same_spin(K, P, Q),
                                 interactions_same_spin(PI<L> - K - Q, PI<L> - P + Q, Q))) {
                    std::cerr << "Interaction parallel is not particle-hole invariant!   " << K << " " << P << " " << Q
                              << ":    " << interactions_same_spin(K, P, Q) << "   "
                              << interactions_same_spin(PI<L> - K - Q, PI<L> - P + Q, Q) << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}

double FlowContainer::band_width() const noexcept {
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for (const auto& point : MomentumGrid::fermion_points()) {
        const double value = dispersion[point.momentum];
        minimum = std::min(minimum, value);
        maximum = std::max(maximum, value);
    }
    return maximum - minimum;
}

double FlowContainer::max_interaction() const noexcept {
    return std::max(interactions_same_spin.norm_inf(), interactions_differing_spin.norm_inf());
}

double FlowContainer::alpha_factor(const momentum_iterator<L>& K,
                                   const momentum_iterator<L>& P,
                                   const momentum_iterator<L>& Q) const noexcept {
#ifdef USE_LW_GENERATOR
        return delta_epsilon(K, P, Q);
#else
        //if (K == -P || Q == PI<L> || Q == Gamma<L>) {
        //    return 0.;
        //}
        //else {
            return sign(delta_epsilon(K, P, Q));
        //}
#endif
}

FlowContainer::coeff_array FlowContainer::epsilon_tilde() const noexcept {
    coeff_array ret;
    for (mom_it K = mom_it::begin(); K != mom_it::end(); ++K) {
        ret[K] = dispersion[K];

        for (const auto& point : MomentumGrid::fermion_points()) {
            const auto P = point.momentum;
            ret[K] -= point.weight * 2. *
                      (interactions_differing_spin(K, P, Gamma<L>) + interactions_same_spin(K, P, Gamma<L>) -
                       interactions_same_spin(K, P, P - K)) *
                      occupation_numbers[P];
        }
    }
    return ret;
};

double FlowContainer::abs_total() const noexcept {
    double val = interactions_same_spin.abs_squared_total() + interactions_differing_spin.abs_squared_total();
    val += dispersion.abs_squared_total();
    return std::sqrt(val);
}

double FlowContainer::norm_inf() const noexcept {
    return std::max(max_interaction(),
                    dispersion.norm_inf());
}

FlowContainer abs(FlowContainer input) noexcept {
    input.interactions_differing_spin.abs_in_place();
    input.interactions_same_spin.abs_in_place();
    input.dispersion.abs_in_place();
    return input;
}

FlowContainer& FlowContainer::operator+=(const FlowContainer& other) noexcept {
    interactions_same_spin += other.interactions_same_spin;
    interactions_differing_spin += other.interactions_differing_spin;
    dispersion += other.dispersion;

    return *this;
}
FlowContainer& FlowContainer::operator-=(const FlowContainer& other) noexcept {
    interactions_same_spin -= other.interactions_same_spin;
    interactions_differing_spin -= other.interactions_differing_spin;
    dispersion -= other.dispersion;

    return *this;
}

FlowContainer& FlowContainer::operator*=(const FlowContainer& other) noexcept {
    interactions_same_spin *= other.interactions_same_spin;
    interactions_differing_spin *= other.interactions_differing_spin;
    dispersion *= other.dispersion;

    return *this;
}
FlowContainer& FlowContainer::operator/=(const FlowContainer& other) noexcept {
    interactions_same_spin /= other.interactions_same_spin;
    interactions_differing_spin /= other.interactions_differing_spin;
    dispersion /= other.dispersion;

    return *this;
}

FlowContainer& FlowContainer::operator*=(const double other) noexcept {
    interactions_same_spin *= other;
    interactions_differing_spin *= other;
    dispersion *= other;

    return *this;
}
FlowContainer& FlowContainer::operator/=(const double other) noexcept {
    interactions_same_spin /= other;
    interactions_differing_spin /= other;
    dispersion /= other;

    return *this;
}

FlowContainer& FlowContainer::operator+=(const double other) noexcept {
    interactions_same_spin += other;
    interactions_differing_spin += other;
    dispersion += other;

    return *this;
}

void to_json(nlohmann::json& j, const FlowContainer& container) noexcept {
    j = nlohmann::json{{"dispersion", container.dispersion.as_array()},
                       {"epsilon_tilde", container.epsilon_tilde()},
                       {"interactions_same_spin", container.interactions_same_spin.as_3D_array()},
                       {"interactions_differing_spin", container.interactions_differing_spin.as_3D_array()}};
}

}  // namespace NickelCUT