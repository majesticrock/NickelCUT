#include "FlowEquation.hpp"

#include "momentum_iterator.hpp"
#include "FlowContainer.hpp"
#include "occupation_numbers.hpp"
#include "../helper_functions.hpp"

#include <omp.h>

namespace NickelCUT::flow {

void FlowEquation::operator()(const FlowContainer& current, FlowContainer& dHdl, const double /*l*/) {
dHdl.reset();
#pragma omp parallel for
for (int K_pos=0; K_pos < N; ++K_pos) {
momentum_iterator<L> K(K_pos);
for (momentum_iterator<L> P = momentum_iterator<L>::begin(); P != momentum_iterator<L>::end(); ++P) {
double nQ_value{};
double one_value{};
for (momentum_iterator<L> Q = momentum_iterator<L>::begin(); Q != momentum_iterator<L>::end(); ++Q) {
one_value += 8.000000 * current.interactions_same_spin(K, P, Q)
	* sign(current.dispersion[-K-Q] + current.dispersion[-P+Q] - current.dispersion[-P] - current.dispersion[-K]) 
	* current.interactions_same_spin(-K-Q, -P+Q, K-P+Q);
one_value -= 8.000000 * current.interactions_same_spin(-Q, -K+P+Q, P+Q)
	* sign(current.dispersion[-K] + current.dispersion[P] - current.dispersion[-Q] - current.dispersion[-K+P+Q]) 
	* current.interactions_same_spin(-K, P, K-Q);
one_value += 8.000000 * current.interactions_differing_spin(-Q, -K+P+Q, P+Q)
	* sign(current.dispersion[-K] + current.dispersion[P] - current.dispersion[-K+P+Q] - current.dispersion[-Q]) 
	* current.interactions_differing_spin(-K, P, P+Q);
one_value -= 8.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.dispersion[-K-Q] + current.dispersion[-P+Q] - current.dispersion[-K] - current.dispersion[-P]) 
	* current.interactions_differing_spin(-K-Q, -P+Q, Q);
one_value += 8.000000 * current.interactions_same_spin(-Q, -K+P+Q, P+Q)
	* sign(current.dispersion[-K] + current.dispersion[P] - current.dispersion[-K+P+Q] - current.dispersion[-Q]) 
	* current.interactions_same_spin(-K, P, P+Q);
one_value -= 8.000000 * current.interactions_same_spin(K, P, Q)
	* sign(current.dispersion[-K-Q] + current.dispersion[-P+Q] - current.dispersion[-K] - current.dispersion[-P]) 
	* current.interactions_same_spin(-K-Q, -P+Q, Q);
nQ_value *= occupation_numbers[Q];
} // Q-loop
dHdl.dispersion[K] += (nQ_value + one_value) * occupation_numbers[P];
} // P-loop
} // K-loop

//--------------------------------------------------------------//
#pragma omp parallel for
for (int K_pos=0; K_pos < N; ++K_pos) {
momentum_iterator<L> K(K_pos);
for (momentum_iterator<L> P = momentum_iterator<L>::begin(); P != momentum_iterator<L>::end(); ++P) {
for (momentum_iterator<L> Q = momentum_iterator<L>::begin(); Q != momentum_iterator<L>::end(); ++Q) {
dHdl.interactions_differing_spin(K, P, Q)+= 4.000000 * sign(current.dispersion[K] + current.dispersion[P] - current.dispersion[K+Q] - current.dispersion[P-Q]) 
	* current.interactions_differing_spin(K, P, Q)
	* current.epsilon_tilde[P-Q];
dHdl.interactions_differing_spin(K, P, Q)-= 4.000000 * sign(current.dispersion[K] + current.dispersion[P] - current.dispersion[K+Q] - current.dispersion[P-Q]) 
	* current.interactions_differing_spin(K, P, Q)
	* current.epsilon_tilde[K];
for (momentum_iterator<L> R = momentum_iterator<L>::begin(); R != momentum_iterator<L>::end(); ++R) {
double nR_value{};
double one_value{};
one_value -= 4.000000 * current.interactions_differing_spin(K, P, R)
	* sign(current.dispersion[K+R] + current.dispersion[P-R] - current.dispersion[K+Q] - current.dispersion[P-Q]) 
	* current.interactions_differing_spin(K+R, P-R, Q-R);
one_value += 4.000000 * current.interactions_differing_spin(-R, K+P+R, K+Q+R)
	* sign(current.dispersion[-K] + current.dispersion[-P] - current.dispersion[R] - current.dispersion[-K-P-R]) 
	* current.interactions_differing_spin(-K, -P, K+R);
dHdl.interactions_differing_spin(K, P, Q) += one_value + occupation_numbers[R] * nR_value;
} // R-loop
} // Q-loop
} // P-loop
} // K-loop

//--------------------------------------------------------------//
#pragma omp parallel for
for (int K_pos=0; K_pos < N; ++K_pos) {
momentum_iterator<L> K(K_pos);
for (momentum_iterator<L> P = momentum_iterator<L>::begin(); P != momentum_iterator<L>::end(); ++P) {
if (K==P) continue; // Pauli principle
for (momentum_iterator<L> Q = momentum_iterator<L>::begin(); Q != momentum_iterator<L>::end(); ++Q) {
if (K+Q==P-Q) continue; // Pauli principle
dHdl.interactions_same_spin(K, P, Q)+= 4.000000 * sign(current.dispersion[K] + current.dispersion[P] - current.dispersion[K+Q] - current.dispersion[P-Q]) 
	* current.interactions_same_spin(K, P, Q)
	* current.epsilon_tilde[P-Q];
dHdl.interactions_same_spin(K, P, Q)-= 4.000000 * sign(current.dispersion[K] + current.dispersion[P] - current.dispersion[K+Q] - current.dispersion[P-Q]) 
	* current.interactions_same_spin(K, P, Q)
	* current.epsilon_tilde[K];
for (momentum_iterator<L> R = momentum_iterator<L>::begin(); R != momentum_iterator<L>::end(); ++R) {
double nR_value{};
double one_value{};
one_value -= 4.000000 * current.interactions_same_spin(K, P, R)
	* sign(current.dispersion[K+R] + current.dispersion[P-R] - current.dispersion[K+Q] - current.dispersion[P-Q]) 
	* current.interactions_same_spin(K+R, P-R, Q-R);
one_value += 4.000000 * current.interactions_same_spin(-R, K+P+R, K+Q+R)
	* sign(current.dispersion[-K] + current.dispersion[-P] - current.dispersion[R] - current.dispersion[-K-P-R]) 
	* current.interactions_same_spin(-K, -P, K+R);
dHdl.interactions_same_spin(K, P, Q) += one_value + occupation_numbers[R] * nR_value;
} // R-loop
} // Q-loop
} // P-loop
} // K-loop
dHdl.interactions_same_spin.symmetrize();
dHdl.interactions_differing_spin.symmetrize();

dHdl.fill_epsilon_tilde();
static int i=0;++i;
if(!(dHdl.is_hermitian() && dHdl.is_inversion_symmetric() && dHdl.is_particle_hole_invariant())) {
	std::cerr << "#" << i << std::endl;
	abort();
}

}
} // namespace NickelCUT::flow