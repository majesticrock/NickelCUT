#include "FlowEquation.hpp"

#include "momentum_iterator.hpp"
#include "FlowContainer.hpp"
#include "occupation_numbers.hpp"
#include "../helper_functions.hpp"

#include <omp.h>

namespace NickelCUT::flow {

void FlowEquation::operator()(const FlowContainer& current, FlowContainer& dHdl, const double /*l*/) {
dHdl.reset();
static InteractionDataFrame alpha_sign_cache;
#pragma omp parallel for
for (int independent_pos = 0; independent_pos < InteractionDataFrame::independent_size(); ++independent_pos) {
const std::size_t tuple_index = InteractionDataFrame::independent_indices()[independent_pos];
momentum_iterator<L> K(static_cast<int>(tuple_index / (N * N)));
momentum_iterator<L> P(static_cast<int>((tuple_index / N) % N));
momentum_iterator<L> Q(static_cast<int>(tuple_index % N));
#ifdef USE_LW_GENERATOR
alpha_sign_cache(K,P,Q) = current.dispersion[K] + current.dispersion[P] - current.dispersion[P-Q] - current.dispersion[K+Q];
#else
alpha_sign_cache(K,P,Q) = sign(current.dispersion[K] + current.dispersion[P] - current.dispersion[P-Q] - current.dispersion[K+Q]);
#endif
}

#pragma omp parallel for
for (int K_pos=0; K_pos < N; ++K_pos) {
momentum_iterator<L> K(K_pos);
for (momentum_iterator<L> P = momentum_iterator<L>::begin(); P != momentum_iterator<L>::end(); ++P) {
double Q_sum_value{};
for (momentum_iterator<L> Q = momentum_iterator<L>::begin(); Q != momentum_iterator<L>::end(); ++Q) {
Q_sum_value -= current.interactions_same_spin(-K, -P, K+Q)
	* alpha_sign_cache(-Q, K+P+Q, K+Q)
	* current.interactions_same_spin(-Q, K+P+Q, K+Q)
	 * hole_numbers[Q] * hole_numbers[K+P+Q];
Q_sum_value += current.interactions_same_spin(-K, -P, K+Q)
	* alpha_sign_cache(-Q, K+P+Q, P+Q)
	* current.interactions_same_spin(-Q, K+P+Q, P+Q)
	 * hole_numbers[Q] * hole_numbers[K+P+Q];
Q_sum_value -= current.interactions_same_spin(-K, K+P+Q, K+P)
	* alpha_sign_cache(-P, -Q, K+P)
	* current.interactions_same_spin(-P, -Q, K+P)
	 * occupation_numbers[Q] * hole_numbers[K+P+Q];
Q_sum_value += current.interactions_same_spin(-K, K+P+Q, K+P)
	* alpha_sign_cache(P, Q, K+Q)
	* current.interactions_same_spin(P, Q, K+Q)
	 * occupation_numbers[Q] * hole_numbers[K+P+Q];
Q_sum_value += current.interactions_same_spin(-P, -Q, K+P)
	* alpha_sign_cache(-K, K+P+Q, K+P)
	* current.interactions_same_spin(-K, K+P+Q, K+P)
	 * occupation_numbers[Q] * hole_numbers[K+P+Q];
Q_sum_value -= current.interactions_same_spin(-P, -Q, K+P)
	* alpha_sign_cache(-K, K+P+Q, K+Q)
	* current.interactions_same_spin(-K, K+P+Q, K+Q)
	 * occupation_numbers[Q] * hole_numbers[K+P+Q];
Q_sum_value += current.interactions_same_spin(-Q, K+P+Q, K+Q)
	* alpha_sign_cache(-K, -P, K+Q)
	* current.interactions_same_spin(-K, -P, K+Q)
	 * hole_numbers[Q] * hole_numbers[K+P+Q];
Q_sum_value -= current.interactions_same_spin(-Q, K+P+Q, K+Q)
	* alpha_sign_cache(-K, -P, -P-Q)
	* current.interactions_same_spin(-K, -P, -P-Q)
	 * hole_numbers[Q] * hole_numbers[K+P+Q];
Q_sum_value -= current.interactions_differing_spin(-K, -P, K+Q)
	* alpha_sign_cache(-Q, K+P+Q, K+Q)
	* current.interactions_differing_spin(-Q, K+P+Q, K+Q)
	 * hole_numbers[Q] * hole_numbers[K+P+Q];
Q_sum_value -= current.interactions_differing_spin(-K, K+P+Q, K+P)
	* alpha_sign_cache(-P, -Q, K+P)
	* current.interactions_differing_spin(-P, -Q, K+P)
	 * occupation_numbers[Q] * hole_numbers[K+P+Q];
Q_sum_value += current.interactions_differing_spin(-P, -Q, K+P)
	* alpha_sign_cache(-K, K+P+Q, K+P)
	* current.interactions_differing_spin(-K, K+P+Q, K+P)
	 * occupation_numbers[Q] * hole_numbers[K+P+Q];
Q_sum_value += current.interactions_differing_spin(-Q, K+P+Q, K+Q)
	* alpha_sign_cache(-K, -P, K+Q)
	* current.interactions_differing_spin(-K, -P, K+Q)
	 * hole_numbers[Q] * hole_numbers[K+P+Q];
} // Q-loop
dHdl.dispersion[K] += 8.000000 * Q_sum_value * occupation_numbers[P];
} // P-loop
} // K-loop

//--------------------------------------------------------------//
#pragma omp parallel for
for (int independent_pos = 0; independent_pos < InteractionDataFrame::independent_size(); ++independent_pos) {
const std::size_t tuple_index = InteractionDataFrame::independent_indices()[independent_pos];
momentum_iterator<L> K(static_cast<int>(tuple_index / (N * N)));
momentum_iterator<L> P(static_cast<int>((tuple_index / N) % N));
momentum_iterator<L> Q(static_cast<int>(tuple_index % N));
dHdl.interactions_differing_spin(K, P, Q) -= alpha_sign_cache(K, P, Q)
	* current.interactions_differing_spin(K, P, Q)
	* current.dispersion[K];
dHdl.interactions_differing_spin(K, P, Q) += alpha_sign_cache(K, P, Q)
	* current.interactions_differing_spin(K, P, Q)
	* current.dispersion[K+Q];
for (momentum_iterator<L> R = momentum_iterator<L>::begin(); R != momentum_iterator<L>::end(); ++R) {
double nR_hQR_value{};
nR_hQR_value -= current.interactions_same_spin(-K, -Q-R, K-R)
	* alpha_sign_cache(-P, -R, Q)
	* current.interactions_differing_spin(-P, -R, Q);
nR_hQR_value += current.interactions_same_spin(K, Q+R, Q)
	* alpha_sign_cache(-P, -R, Q)
	* current.interactions_differing_spin(-P, -R, Q);
nR_hQR_value -= current.interactions_same_spin(-P, -R, Q)
	* alpha_sign_cache(K, Q+R, Q)
	* current.interactions_differing_spin(K, Q+R, Q);
nR_hQR_value += current.interactions_same_spin(-P, -R, P-Q-R)
	* alpha_sign_cache(K, Q+R, Q)
	* current.interactions_differing_spin(K, Q+R, Q);
dHdl.interactions_differing_spin(K, P, Q) += current.interactions_differing_spin(-K, -P, K+R)
	* alpha_sign_cache(-R, K+P+R, K+Q+R)
	* current.interactions_differing_spin(-R, K+P+R, K+Q+R)
	 * occupation_numbers[R] * occupation_numbers[K+P+R];
dHdl.interactions_differing_spin(K, P, Q) -= current.interactions_differing_spin(-K, -P, K+R)
	* alpha_sign_cache(-R, K+P+R, K+Q+R)
	* current.interactions_differing_spin(-R, K+P+R, K+Q+R)
	 * hole_numbers[R] * hole_numbers[K+P+R];
dHdl.interactions_differing_spin(K, P, Q) += 2.000000 * current.interactions_differing_spin(-K, -R, P-Q-R)
	* alpha_sign_cache(-P, -K+P-Q-R, P-R)
	* current.interactions_differing_spin(-P, -K+P-Q-R, P-R)
	 * occupation_numbers[R] * hole_numbers[K-P+Q+R];
nR_hQR_value += current.interactions_differing_spin(K, Q+R, Q)
	* alpha_sign_cache(-P, -R, Q)
	* current.interactions_same_spin(-P, -R, Q);
nR_hQR_value -= current.interactions_differing_spin(K, Q+R, Q)
	* alpha_sign_cache(-P, -R, P-Q-R)
	* current.interactions_same_spin(-P, -R, P-Q-R);
nR_hQR_value += current.interactions_differing_spin(-P, -R, Q)
	* alpha_sign_cache(-K, -Q-R, K-R)
	* current.interactions_same_spin(-K, -Q-R, K-R);
nR_hQR_value -= current.interactions_differing_spin(-P, -R, Q)
	* alpha_sign_cache(K, Q+R, Q)
	* current.interactions_same_spin(K, Q+R, Q);
dHdl.interactions_differing_spin(K, P, Q) -= 2.000000 * current.interactions_differing_spin(-P, -K+P-Q-R, P-R)
	* alpha_sign_cache(-K, -R, P-Q-R)
	* current.interactions_differing_spin(-K, -R, P-Q-R)
	 * occupation_numbers[R] * hole_numbers[K-P+Q+R];
dHdl.interactions_differing_spin(K, P, Q) -= current.interactions_differing_spin(-R, K+P+R, K+Q+R)
	* alpha_sign_cache(-K, -P, K+R)
	* current.interactions_differing_spin(-K, -P, K+R)
	 * occupation_numbers[R] * occupation_numbers[K+P+R];
dHdl.interactions_differing_spin(K, P, Q) += current.interactions_differing_spin(-R, K+P+R, K+Q+R)
	* alpha_sign_cache(-K, -P, K+R)
	* current.interactions_differing_spin(-K, -P, K+R)
	 * hole_numbers[R] * hole_numbers[K+P+R];
dHdl.interactions_differing_spin(K, P, Q) += 2. * nR_hQR_value * occupation_numbers[R] * hole_numbers[Q+R];
} // R-loop
dHdl.interactions_differing_spin(K, P, Q) *= 4.000000;
}

//--------------------------------------------------------------//
#pragma omp parallel for
for (int independent_pos = 0; independent_pos < InteractionDataFrame::independent_size(); ++independent_pos) {
const std::size_t tuple_index = InteractionDataFrame::independent_indices()[independent_pos];
momentum_iterator<L> K(static_cast<int>(tuple_index / (N * N)));
momentum_iterator<L> P(static_cast<int>((tuple_index / N) % N));
momentum_iterator<L> Q(static_cast<int>(tuple_index % N));
if (K==P) continue; // Pauli principle
if (K+Q==P-Q) continue; // Pauli principle
dHdl.interactions_same_spin(K, P, Q) -= alpha_sign_cache(K, P, Q)
	* current.interactions_same_spin(K, P, Q)
	* current.dispersion[K];
dHdl.interactions_same_spin(K, P, Q) += alpha_sign_cache(K, P, Q)
	* current.interactions_same_spin(K, P, Q)
	* current.dispersion[K+Q];
for (momentum_iterator<L> R = momentum_iterator<L>::begin(); R != momentum_iterator<L>::end(); ++R) {
double nR_hQR_value{};
dHdl.interactions_same_spin(K, P, Q) += current.interactions_same_spin(-K, -P, K+R)
	* alpha_sign_cache(-R, K+P+R, K+Q+R)
	* current.interactions_same_spin(-R, K+P+R, K+Q+R)
	 * occupation_numbers[R] * occupation_numbers[K+P+R];
dHdl.interactions_same_spin(K, P, Q) -= current.interactions_same_spin(-K, -P, K+R)
	* alpha_sign_cache(-R, K+P+R, K+Q+R)
	* current.interactions_same_spin(-R, K+P+R, K+Q+R)
	 * hole_numbers[R] * hole_numbers[K+P+R];
dHdl.interactions_same_spin(K, P, Q) += 2.000000 * current.interactions_same_spin(-K, -R, P-Q-R)
	* alpha_sign_cache(-P, -K+P-Q-R, P-R)
	* current.interactions_same_spin(-P, -K+P-Q-R, P-R)
	 * occupation_numbers[R] * hole_numbers[K-P+Q+R];
nR_hQR_value -= current.interactions_same_spin(-K, -Q-R, K-R)
	* alpha_sign_cache(-P, -R, Q)
	* current.interactions_same_spin(-P, -R, Q);
nR_hQR_value += current.interactions_same_spin(K, Q+R, Q)
	* alpha_sign_cache(-P, -R, Q)
	* current.interactions_same_spin(-P, -R, Q);
nR_hQR_value -= current.interactions_same_spin(K, Q+R, Q)
	* alpha_sign_cache(-P, -R, P-Q-R)
	* current.interactions_same_spin(-P, -R, P-Q-R);
nR_hQR_value += current.interactions_same_spin(-P, -R, Q)
	* alpha_sign_cache(-K, -Q-R, K-R)
	* current.interactions_same_spin(-K, -Q-R, K-R);
nR_hQR_value -= current.interactions_same_spin(-P, -R, Q)
	* alpha_sign_cache(K, Q+R, Q)
	* current.interactions_same_spin(K, Q+R, Q);
nR_hQR_value += current.interactions_same_spin(-P, -R, P-Q-R)
	* alpha_sign_cache(K, Q+R, Q)
	* current.interactions_same_spin(K, Q+R, Q);
dHdl.interactions_same_spin(K, P, Q) -= 2.000000 * current.interactions_same_spin(-P, -K+P-Q-R, P-R)
	* alpha_sign_cache(-K, -R, P-Q-R)
	* current.interactions_same_spin(-K, -R, P-Q-R)
	 * occupation_numbers[R] * hole_numbers[K-P+Q+R];
dHdl.interactions_same_spin(K, P, Q) -= current.interactions_same_spin(-R, K+P+R, K+Q+R)
	* alpha_sign_cache(-K, -P, K+R)
	* current.interactions_same_spin(-K, -P, K+R)
	 * occupation_numbers[R] * occupation_numbers[K+P+R];
dHdl.interactions_same_spin(K, P, Q) += current.interactions_same_spin(-R, K+P+R, K+Q+R)
	* alpha_sign_cache(-K, -P, K+R)
	* current.interactions_same_spin(-K, -P, K+R)
	 * hole_numbers[R] * hole_numbers[K+P+R];
nR_hQR_value += current.interactions_differing_spin(K, Q+R, Q)
	* alpha_sign_cache(-P, -R, Q)
	* current.interactions_differing_spin(-P, -R, Q);
nR_hQR_value -= current.interactions_differing_spin(-P, -R, Q)
	* alpha_sign_cache(K, Q+R, Q)
	* current.interactions_differing_spin(K, Q+R, Q);
dHdl.interactions_same_spin(K, P, Q) += 2. * nR_hQR_value * occupation_numbers[R] * hole_numbers[Q+R];
} // R-loop
dHdl.interactions_same_spin(K, P, Q) *= 4.000000;
}
dHdl.interactions_same_spin.antisymmetrize();
dHdl.interactions_same_spin.symmetrize();
dHdl.interactions_differing_spin.symmetrize();
dHdl.interactions_same_spin.clear_noise();
dHdl.interactions_differing_spin.clear_noise();
}
} // namespace NickelCUT::flow