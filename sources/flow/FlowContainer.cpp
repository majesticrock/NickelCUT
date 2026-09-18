#include "FlowContainer.hpp"
#include "Model.hpp"
#include "momentum_iterator.hpp"
#include "occupation_numbers.hpp"
#include "../helper_functions.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <iostream>

namespace NickelCUT::flow {

FlowContainer::FlowContainer() noexcept
    : interactions_same_spin(), interactions_differing_spin()
{
    dispersion.fill(double{});
    epsilon_tilde.fill(double{});
}

FlowContainer::FlowContainer(const Model& model) noexcept 
    : interactions_same_spin(), 
    // Factor 1/2 to account for the extra spin summation compared to the standard Hubbard model
    interactions_differing_spin(0.5 * model.U_0 / N)
{
    for (mom_it K = mom_it::begin(); K != mom_it::end(); ++K) {
        dispersion[K] = model.epsilon_0(K.get_kx(), K.get_ky());// + 0.25 * model.U_0;
    }
    fill_epsilon_tilde();
}

void FlowContainer::print_memory_allocation() const
{
    const std::size_t allocated_bytes = sizeof(FlowContainer)
        + interactions_same_spin.get_data().capacity() * sizeof(double)
        + interactions_differing_spin.get_data().capacity() * sizeof(double);

    std::cout << "FlowContainer constructed, allocated memory: "
              << allocated_bytes / (1024) << " kb" << std::endl;
}

void FlowContainer::fill(double value) noexcept
{
    interactions_same_spin.fill(value);
    interactions_differing_spin.fill(value);
    dispersion.fill(value);
    epsilon_tilde.fill(value);
}

void FlowContainer::reset() noexcept
{
    fill(double{});
}

bool FlowContainer::contains_nan_or_inf() const noexcept {
    for (const auto& val : interactions_same_spin.get_data()) {
        if (std::isnan(val) || std::isinf(val)) return true;
    }
    for (const auto& val : interactions_differing_spin.get_data()) {
        if (std::isnan(val) || std::isinf(val)) return true;
    }
    for (const auto& val : dispersion) {
        if (std::isnan(val) || std::isinf(val)) return true;
    }
    for (const auto& val : epsilon_tilde) {
        if (std::isnan(val) || std::isinf(val)) return true;
    }
    return false;
}

double FlowContainer::residual_offdiagonality() const noexcept {
    return std::sqrt(interactions_same_spin.abs_squared_total() + interactions_differing_spin.abs_squared_total());
}

bool FlowContainer::is_inversion_symmetric() const noexcept
{
    for (mom_it K = mom_it::begin(); K != mom_it::end(); ++K) {
        if (!float_equal(dispersion[K], dispersion[-K])) {
            std::cerr << "Dispersion is not inversion symmetric!" << std::endl;
            return false;
        }
    }

    for (mom_it K = mom_it::begin(); K != mom_it::end(); ++K) {
        for (mom_it P = mom_it::begin(); P != mom_it::end(); ++P) {
            for (mom_it Q = mom_it::begin(); Q != mom_it::end(); ++Q) {
                if(!float_equal(interactions_differing_spin(K, P, Q), interactions_differing_spin(-K, -P, -Q))) {
                    std::cerr << "Interaction antiparallel is not inversion symmetric!" << std::endl;
                    return false;
                }
            }
        }
    }

    for (mom_it K = mom_it::begin(); K != mom_it::end(); ++K) {
        for (mom_it P = mom_it::begin(); P != mom_it::end(); ++P) {
            for (mom_it Q = mom_it::begin(); Q != mom_it::end(); ++Q) {
                if(!float_equal(interactions_same_spin(K, P, Q), interactions_same_spin(-K, -P, -Q))) {
                    std::cerr << "Interaction parallel is not inversion symmetric!" << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}

bool FlowContainer::is_hermitian() const noexcept
{
    for (mom_it K = mom_it::begin(); K != mom_it::end(); ++K) {
        for (mom_it P = mom_it::begin(); P != mom_it::end(); ++P) {
            for (mom_it Q = mom_it::begin(); Q != mom_it::end(); ++Q) {
                if(!float_equal(interactions_differing_spin(K, P, Q), interactions_differing_spin(K+Q, P-Q, -Q))) {
                    std::cerr << "Interaction antiparallel is not Hermitian!" << std::endl;
                    return false;
                }
            }
        }
    }

    for (mom_it K = mom_it::begin(); K != mom_it::end(); ++K) {
        for (mom_it P = mom_it::begin(); P != mom_it::end(); ++P) {
            for (mom_it Q = mom_it::begin(); Q != mom_it::end(); ++Q) {
                if(!float_equal(interactions_same_spin(K, P, Q), interactions_same_spin(K+Q, P-Q, -Q))) {
                    std::cerr << "Interaction parallel is not Hermitian!" << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}

bool FlowContainer::is_particle_hole_invariant() const noexcept
{
    for (mom_it K = mom_it::begin(); K != mom_it::end(); ++K) {
        double interaction_contribution = 0.;
        for (mom_it P = mom_it::begin(); P != mom_it::end(); ++P) {
            interaction_contribution += interactions_same_spin(PI<L>-K, PI<L>-K+P, P);
            interaction_contribution -= interactions_same_spin(PI<L>-K, P, 0);
            interaction_contribution -= interactions_differing_spin(PI<L>-K, P, 0);
        }
        interaction_contribution *= 2.;
        if (!float_equal(epsilon_tilde[K] + epsilon_tilde[PI<L>-K], interaction_contribution)) {
            std::cerr << "Dispersion is not particle-hole invariant!       "  << K << ":\t"
                << epsilon_tilde[K] << "    " << epsilon_tilde[PI<L>-K]
                << "    " << interaction_contribution
                << "    " << interaction_contribution-epsilon_tilde[PI<L>-K]-epsilon_tilde[K] << std::endl;
            //return false;
        }
    }

    for (mom_it K = mom_it::begin(); K != mom_it::end(); ++K) {
        for (mom_it P = mom_it::begin(); P != mom_it::end(); ++P) {
            for (mom_it Q = mom_it::begin(); Q != mom_it::end(); ++Q) {
                if(!float_equal(interactions_differing_spin(K, P, Q), interactions_differing_spin(PI<L>-K-Q, PI<L>-P+Q, Q))) {
                    std::cerr << "Interaction antiparallel is not particle-hole invariant!   " 
                        << K << " " << P << " " << Q << ":    "
                        << interactions_differing_spin(K, P, Q) << "   " 
                        << interactions_differing_spin(PI<L>-K-Q, PI<L>-P+Q, Q) << std::endl;

                    return false;
                }
            }
        }
    }

    for (mom_it K = mom_it::begin(); K != mom_it::end(); ++K) {
        for (mom_it P = mom_it::begin(); P != mom_it::end(); ++P) {
            for (mom_it Q = mom_it::begin(); Q != mom_it::end(); ++Q) {
                if(!float_equal(interactions_same_spin(K, P, Q), interactions_same_spin(PI<L>-K-Q, PI<L>-P+Q, Q))) {
                    std::cerr << "Interaction parallel is not particle-hole invariant!   " 
                        << K << " " << P << " " << Q << ":    "
                        << interactions_same_spin(K, P, Q) << "   " 
                        << interactions_same_spin(PI<L>-K-Q, PI<L>-P+Q, Q) << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}

void FlowContainer::fill_epsilon_tilde() noexcept {
    for (mom_it K = mom_it::begin(); K != mom_it::end(); ++K) {
        epsilon_tilde[K] = dispersion[K];

        for (mom_it P = mom_it::begin(); P != mom_it::end(); ++P) {
            epsilon_tilde[K] -= 2. * (interactions_differing_spin(K, P, Gamma<L>)
                                    + interactions_same_spin(K, P, Gamma<L>)
                                    - interactions_same_spin(K, P, P-K)
                                 ) * occupation_numbers[P];
        }
    }
};

double FlowContainer::abs_total() const noexcept
{
    double val = interactions_same_spin.abs_squared_total() + interactions_differing_spin.abs_squared_total();
    for (const auto& element : dispersion) {
        val += element*element;
    }
    return std::sqrt(val);
}

double FlowContainer::norm_inf() const noexcept
{
    double val = interactions_same_spin.norm_inf();

    double comp = interactions_differing_spin.norm_inf();
    if(val < comp) val = comp;

    comp = std::abs(*std::max_element(dispersion.begin(), dispersion.end(), LessThanAbs()));
    if(val < comp) val = comp;

    return val;
}

FlowContainer abs(FlowContainer input) noexcept
{
    input.interactions_differing_spin.abs_in_place();
    input.interactions_same_spin.abs_in_place();
    for (auto& val : input.dispersion) {
        val = std::abs(val);
    }
    return input;
}

FlowContainer& FlowContainer::operator+=(const FlowContainer& other) noexcept
{
    interactions_same_spin += other.interactions_same_spin;
    interactions_differing_spin += other.interactions_differing_spin;
    for (int i=0; i<N; ++i) {
        dispersion[i] += other.dispersion[i];
    }

    fill_epsilon_tilde();

    return *this;
}
FlowContainer& FlowContainer::operator-=(const FlowContainer& other) noexcept
{
    interactions_same_spin -= other.interactions_same_spin;
    interactions_differing_spin -= other.interactions_differing_spin;
    for (int i=0; i<N; ++i) {
        dispersion[i] -= other.dispersion[i];
    }

    fill_epsilon_tilde();

    return *this;
}

FlowContainer& FlowContainer::operator*=(const FlowContainer& other) noexcept
{
    interactions_same_spin *= other.interactions_same_spin;
    interactions_differing_spin *= other.interactions_differing_spin;
    for (int i=0; i<N; ++i) {
        dispersion[i] *= other.dispersion[i];
    }

    fill_epsilon_tilde();

    return *this;
}
FlowContainer& FlowContainer::operator/=(const FlowContainer& other) noexcept
{
    interactions_same_spin /= other.interactions_same_spin;
    interactions_differing_spin /= other.interactions_differing_spin;
    for (int i=0; i<N; ++i) {
        dispersion[i] /= other.dispersion[i];
    }

    fill_epsilon_tilde();

    return *this;
}

FlowContainer& FlowContainer::operator*=(const double other) noexcept
{
    interactions_same_spin *= other;
    interactions_differing_spin *= other;
    for (int i=0; i<N; ++i) {
        dispersion[i] *= other;
    }

    fill_epsilon_tilde();

    return *this;
}
FlowContainer& FlowContainer::operator/=(const double other) noexcept
{
    interactions_same_spin /= other;
    interactions_differing_spin /= other;
    for (int i=0; i<N; ++i) {
        dispersion[i] /= other;
    }

    fill_epsilon_tilde();

    return *this;
}

FlowContainer& FlowContainer::operator+=(const double other) noexcept
{
    interactions_same_spin += other;
    interactions_differing_spin += other;
    for (int i=0; i<N; ++i) {
        dispersion[i] += other;
    }

    fill_epsilon_tilde();

    return *this;
}

void to_json(nlohmann::json& j, const FlowContainer& container) noexcept
{
    j = nlohmann::json{
        { "dispersion",                  container.dispersion                                },
        { "epsilon_tilde",               container.epsilon_tilde                             },
        { "interactions_same_spin",      container.interactions_same_spin.as_3D_array()      },
        { "interactions_differing_spin", container.interactions_differing_spin.as_3D_array() }
    };
}

} // namespace NickelCUT::flow