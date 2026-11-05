#include "FlowEquation.hpp"

#include "../helper_functions.hpp"
#include "../shared/FlowContainer.hpp"
#include "../shared/MomentumGrid.hpp"
#include "../momentum_iterator.hpp"
#include "../shared/occupation_numbers.hpp"

#include <omp.h>

namespace NickelCUT::flow {

void FlowEquation::operator()(const FlowContainer& current, FlowContainer& dHdl, const double /*l*/) {
    dHdl.reset();
    static InteractionDataFrame alpha_sign_cache;
    const auto& active_interaction_tuples = InteractionDataFrame::active_independent_indices();
    const auto& active_fermion_indices = SingleParticleDataFrame::active_independent_indices();
#pragma omp parallel for
    for (int independent_pos = 0;
         independent_pos < static_cast<int>(active_interaction_tuples.size());
         ++independent_pos) {
        const std::size_t tuple_index = active_interaction_tuples[independent_pos];
        momentum_iterator<L> K(static_cast<int>(tuple_index / (N * N)));
        momentum_iterator<L> P(static_cast<int>((tuple_index / N) % N));
        momentum_iterator<L> Q(static_cast<int>(tuple_index % N));
        alpha_sign_cache(K, P, Q) = current.alpha_factor(K, P, Q);
    }

#pragma omp parallel for
    for (int K_pos = 0; K_pos < static_cast<int>(active_fermion_indices.size()); ++K_pos) {
        const std::size_t k_index = active_fermion_indices[K_pos];
        momentum_iterator<L> K(static_cast<int>(k_index));
        for (const auto& P_point : MomentumGrid::fermion_points()) {
            const auto P = P_point.momentum;
            double Q_sum_value{};
            for (const auto& Q_point : MomentumGrid::transfer_points()) {
                const auto Q = Q_point.momentum;
                const double q_sum_before = Q_sum_value;
                Q_sum_value -= current.interactions_same_spin(-K, -P, K + Q) * alpha_sign_cache(-Q, K + P + Q, K + Q) *
                               current.interactions_same_spin(-Q, K + P + Q, K + Q) * hole_numbers[Q] *
                               hole_numbers[K + P + Q];
                Q_sum_value += current.interactions_same_spin(-K, -P, K + Q) * alpha_sign_cache(-Q, K + P + Q, P + Q) *
                               current.interactions_same_spin(-Q, K + P + Q, P + Q) * hole_numbers[Q] *
                               hole_numbers[K + P + Q];
                Q_sum_value -= current.interactions_same_spin(-K, K + P + Q, K + P) * alpha_sign_cache(-P, -Q, K + P) *
                               current.interactions_same_spin(-P, -Q, K + P) * occupation_numbers[Q] *
                               hole_numbers[K + P + Q];
                Q_sum_value += current.interactions_same_spin(-K, K + P + Q, K + P) * alpha_sign_cache(P, Q, K + Q) *
                               current.interactions_same_spin(P, Q, K + Q) * occupation_numbers[Q] *
                               hole_numbers[K + P + Q];
                Q_sum_value += current.interactions_same_spin(-P, -Q, K + P) * alpha_sign_cache(-K, K + P + Q, K + P) *
                               current.interactions_same_spin(-K, K + P + Q, K + P) * occupation_numbers[Q] *
                               hole_numbers[K + P + Q];
                Q_sum_value -= current.interactions_same_spin(-P, -Q, K + P) * alpha_sign_cache(-K, K + P + Q, K + Q) *
                               current.interactions_same_spin(-K, K + P + Q, K + Q) * occupation_numbers[Q] *
                               hole_numbers[K + P + Q];
                Q_sum_value += current.interactions_same_spin(-Q, K + P + Q, K + Q) * alpha_sign_cache(-K, -P, K + Q) *
                               current.interactions_same_spin(-K, -P, K + Q) * hole_numbers[Q] *
                               hole_numbers[K + P + Q];
                Q_sum_value -= current.interactions_same_spin(-Q, K + P + Q, K + Q) * alpha_sign_cache(-K, -P, -P - Q) *
                               current.interactions_same_spin(-K, -P, -P - Q) * hole_numbers[Q] *
                               hole_numbers[K + P + Q];
                Q_sum_value -= current.interactions_differing_spin(-K, -P, K + Q) *
                               alpha_sign_cache(-Q, K + P + Q, K + Q) *
                               current.interactions_differing_spin(-Q, K + P + Q, K + Q) * hole_numbers[Q] *
                               hole_numbers[K + P + Q];
                Q_sum_value -= current.interactions_differing_spin(-K, K + P + Q, K + P) *
                               alpha_sign_cache(-P, -Q, K + P) * current.interactions_differing_spin(-P, -Q, K + P) *
                               occupation_numbers[Q] * hole_numbers[K + P + Q];
                Q_sum_value += current.interactions_differing_spin(-P, -Q, K + P) *
                               alpha_sign_cache(-K, K + P + Q, K + P) *
                               current.interactions_differing_spin(-K, K + P + Q, K + P) * occupation_numbers[Q] *
                               hole_numbers[K + P + Q];
                Q_sum_value += current.interactions_differing_spin(-Q, K + P + Q, K + Q) *
                               alpha_sign_cache(-K, -P, K + Q) * current.interactions_differing_spin(-K, -P, K + Q) *
                               hole_numbers[Q] * hole_numbers[K + P + Q];
                Q_sum_value = q_sum_before + (Q_sum_value - q_sum_before) * Q_point.weight;
            }  // Q-loop
            dHdl.dispersion[K] += 8.000000 * Q_sum_value * occupation_numbers[P] * P_point.weight;
        }  // P-loop
    }  // K-loop

//--------------------------------------------------------------//
#pragma omp parallel for
    for (int independent_pos = 0;
         independent_pos < static_cast<int>(active_interaction_tuples.size());
         ++independent_pos) {
        const std::size_t tuple_index = active_interaction_tuples[independent_pos];
        momentum_iterator<L> K(static_cast<int>(tuple_index / (N * N)));
        momentum_iterator<L> P(static_cast<int>((tuple_index / N) % N));
        momentum_iterator<L> Q(static_cast<int>(tuple_index % N));
        dHdl.interactions_differing_spin(K, P, Q) -=
            alpha_sign_cache(K, P, Q) * current.interactions_differing_spin(K, P, Q) * current.dispersion[K];
        dHdl.interactions_differing_spin(K, P, Q) +=
            alpha_sign_cache(K, P, Q) * current.interactions_differing_spin(K, P, Q) * current.dispersion[K + Q];
        for (const auto& R_point : MomentumGrid::fermion_points()) {
            const auto R = R_point.momentum;
            const double r_sum_before = dHdl.interactions_differing_spin(K, P, Q);
            double nR_hQR_value{};
            nR_hQR_value -= current.interactions_same_spin(-K, -Q - R, K - R) * alpha_sign_cache(-P, -R, Q) *
                            current.interactions_differing_spin(-P, -R, Q);
            nR_hQR_value += current.interactions_same_spin(K, Q + R, Q) * alpha_sign_cache(-P, -R, Q) *
                            current.interactions_differing_spin(-P, -R, Q);
            nR_hQR_value -= current.interactions_same_spin(-P, -R, Q) * alpha_sign_cache(K, Q + R, Q) *
                            current.interactions_differing_spin(K, Q + R, Q);
            nR_hQR_value += current.interactions_same_spin(-P, -R, P - Q - R) * alpha_sign_cache(K, Q + R, Q) *
                            current.interactions_differing_spin(K, Q + R, Q);
            dHdl.interactions_differing_spin(K, P, Q) += current.interactions_differing_spin(-K, -P, K + R) *
                                                         alpha_sign_cache(-R, K + P + R, K + Q + R) *
                                                         current.interactions_differing_spin(-R, K + P + R, K + Q + R) *
                                                         occupation_numbers[R] * occupation_numbers[K + P + R];
            dHdl.interactions_differing_spin(K, P, Q) -= current.interactions_differing_spin(-K, -P, K + R) *
                                                         alpha_sign_cache(-R, K + P + R, K + Q + R) *
                                                         current.interactions_differing_spin(-R, K + P + R, K + Q + R) *
                                                         hole_numbers[R] * hole_numbers[K + P + R];
            dHdl.interactions_differing_spin(K, P, Q) +=
                2.000000 * current.interactions_differing_spin(-K, -R, P - Q - R) *
                alpha_sign_cache(-P, -K + P - Q - R, P - R) *
                current.interactions_differing_spin(-P, -K + P - Q - R, P - R) * occupation_numbers[R] *
                hole_numbers[K - P + Q + R];
            nR_hQR_value += current.interactions_differing_spin(K, Q + R, Q) * alpha_sign_cache(-P, -R, Q) *
                            current.interactions_same_spin(-P, -R, Q);
            nR_hQR_value -= current.interactions_differing_spin(K, Q + R, Q) * alpha_sign_cache(-P, -R, P - Q - R) *
                            current.interactions_same_spin(-P, -R, P - Q - R);
            nR_hQR_value += current.interactions_differing_spin(-P, -R, Q) * alpha_sign_cache(-K, -Q - R, K - R) *
                            current.interactions_same_spin(-K, -Q - R, K - R);
            nR_hQR_value -= current.interactions_differing_spin(-P, -R, Q) * alpha_sign_cache(K, Q + R, Q) *
                            current.interactions_same_spin(K, Q + R, Q);
            dHdl.interactions_differing_spin(K, P, Q) -=
                2.000000 * current.interactions_differing_spin(-P, -K + P - Q - R, P - R) *
                alpha_sign_cache(-K, -R, P - Q - R) * current.interactions_differing_spin(-K, -R, P - Q - R) *
                occupation_numbers[R] * hole_numbers[K - P + Q + R];
            dHdl.interactions_differing_spin(K, P, Q) -= current.interactions_differing_spin(-R, K + P + R, K + Q + R) *
                                                         alpha_sign_cache(-K, -P, K + R) *
                                                         current.interactions_differing_spin(-K, -P, K + R) *
                                                         occupation_numbers[R] * occupation_numbers[K + P + R];
            dHdl.interactions_differing_spin(K, P, Q) +=
                current.interactions_differing_spin(-R, K + P + R, K + Q + R) * alpha_sign_cache(-K, -P, K + R) *
                current.interactions_differing_spin(-K, -P, K + R) * hole_numbers[R] * hole_numbers[K + P + R];
            dHdl.interactions_differing_spin(K, P, Q) +=
                2. * nR_hQR_value * occupation_numbers[R] * hole_numbers[Q + R];
            const double r_sum_after = dHdl.interactions_differing_spin(K, P, Q);
            dHdl.interactions_differing_spin(K, P, Q) =
                r_sum_before + (r_sum_after - r_sum_before) * R_point.weight;
        }  // R-loop
        dHdl.interactions_differing_spin(K, P, Q) *= 4.000000;
    }

//--------------------------------------------------------------//
#pragma omp parallel for
    for (int independent_pos = 0;
         independent_pos < static_cast<int>(active_interaction_tuples.size());
         ++independent_pos) {
        const std::size_t tuple_index = active_interaction_tuples[independent_pos];
        momentum_iterator<L> K(static_cast<int>(tuple_index / (N * N)));
        momentum_iterator<L> P(static_cast<int>((tuple_index / N) % N));
        momentum_iterator<L> Q(static_cast<int>(tuple_index % N));
        if (K == P)
            continue;  // Pauli principle
        if (K + Q == P - Q)
            continue;  // Pauli principle
        dHdl.interactions_same_spin(K, P, Q) -=
            alpha_sign_cache(K, P, Q) * current.interactions_same_spin(K, P, Q) * current.dispersion[K];
        dHdl.interactions_same_spin(K, P, Q) +=
            alpha_sign_cache(K, P, Q) * current.interactions_same_spin(K, P, Q) * current.dispersion[K + Q];
        for (const auto& R_point : MomentumGrid::fermion_points()) {
            const auto R = R_point.momentum;
            const double r_sum_before = dHdl.interactions_same_spin(K, P, Q);
            double nR_hQR_value{};
            dHdl.interactions_same_spin(K, P, Q) += current.interactions_same_spin(-K, -P, K + R) *
                                                    alpha_sign_cache(-R, K + P + R, K + Q + R) *
                                                    current.interactions_same_spin(-R, K + P + R, K + Q + R) *
                                                    occupation_numbers[R] * occupation_numbers[K + P + R];
            dHdl.interactions_same_spin(K, P, Q) -=
                current.interactions_same_spin(-K, -P, K + R) * alpha_sign_cache(-R, K + P + R, K + Q + R) *
                current.interactions_same_spin(-R, K + P + R, K + Q + R) * hole_numbers[R] * hole_numbers[K + P + R];
            dHdl.interactions_same_spin(K, P, Q) += 2.000000 * current.interactions_same_spin(-K, -R, P - Q - R) *
                                                    alpha_sign_cache(-P, -K + P - Q - R, P - R) *
                                                    current.interactions_same_spin(-P, -K + P - Q - R, P - R) *
                                                    occupation_numbers[R] * hole_numbers[K - P + Q + R];
            nR_hQR_value -= current.interactions_same_spin(-K, -Q - R, K - R) * alpha_sign_cache(-P, -R, Q) *
                            current.interactions_same_spin(-P, -R, Q);
            nR_hQR_value += current.interactions_same_spin(K, Q + R, Q) * alpha_sign_cache(-P, -R, Q) *
                            current.interactions_same_spin(-P, -R, Q);
            nR_hQR_value -= current.interactions_same_spin(K, Q + R, Q) * alpha_sign_cache(-P, -R, P - Q - R) *
                            current.interactions_same_spin(-P, -R, P - Q - R);
            nR_hQR_value += current.interactions_same_spin(-P, -R, Q) * alpha_sign_cache(-K, -Q - R, K - R) *
                            current.interactions_same_spin(-K, -Q - R, K - R);
            nR_hQR_value -= current.interactions_same_spin(-P, -R, Q) * alpha_sign_cache(K, Q + R, Q) *
                            current.interactions_same_spin(K, Q + R, Q);
            nR_hQR_value += current.interactions_same_spin(-P, -R, P - Q - R) * alpha_sign_cache(K, Q + R, Q) *
                            current.interactions_same_spin(K, Q + R, Q);
            dHdl.interactions_same_spin(K, P, Q) -=
                2.000000 * current.interactions_same_spin(-P, -K + P - Q - R, P - R) *
                alpha_sign_cache(-K, -R, P - Q - R) * current.interactions_same_spin(-K, -R, P - Q - R) *
                occupation_numbers[R] * hole_numbers[K - P + Q + R];
            dHdl.interactions_same_spin(K, P, Q) -=
                current.interactions_same_spin(-R, K + P + R, K + Q + R) * alpha_sign_cache(-K, -P, K + R) *
                current.interactions_same_spin(-K, -P, K + R) * occupation_numbers[R] * occupation_numbers[K + P + R];
            dHdl.interactions_same_spin(K, P, Q) +=
                current.interactions_same_spin(-R, K + P + R, K + Q + R) * alpha_sign_cache(-K, -P, K + R) *
                current.interactions_same_spin(-K, -P, K + R) * hole_numbers[R] * hole_numbers[K + P + R];
            nR_hQR_value += current.interactions_differing_spin(K, Q + R, Q) * alpha_sign_cache(-P, -R, Q) *
                            current.interactions_differing_spin(-P, -R, Q);
            nR_hQR_value -= current.interactions_differing_spin(-P, -R, Q) * alpha_sign_cache(K, Q + R, Q) *
                            current.interactions_differing_spin(K, Q + R, Q);
            dHdl.interactions_same_spin(K, P, Q) += 2. * nR_hQR_value * occupation_numbers[R] * hole_numbers[Q + R];
            const double r_sum_after = dHdl.interactions_same_spin(K, P, Q);
            dHdl.interactions_same_spin(K, P, Q) =
                r_sum_before + (r_sum_after - r_sum_before) * R_point.weight;
        }  // R-loop
        dHdl.interactions_same_spin(K, P, Q) *= 4.000000;
    }
    dHdl.interactions_same_spin.antisymmetrize();
    dHdl.interactions_same_spin.symmetrize();
    dHdl.interactions_differing_spin.symmetrize();
    dHdl.interactions_same_spin.clear_noise();
    dHdl.interactions_differing_spin.clear_noise();
}
}  // namespace NickelCUT::flow