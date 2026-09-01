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
one_value -= 8.000000 * current.interactions_same_spin(-Q, -K+P+Q, P+Q)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[P] - current.epsilon_tilde[-Q] - current.epsilon_tilde[-K+P+Q]) 
	* current.interactions_same_spin(-K, P, K-Q);
one_value += 8.000000 * current.interactions_same_spin(K, P, Q)
	* sign(current.epsilon_tilde[-K-Q] + current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-P] - current.epsilon_tilde[-K]) 
	* current.interactions_same_spin(-K-Q, -P+Q, K-P+Q);
one_value += 8.000000 * current.interactions_differing_spin(-Q, -K+P+Q, P+Q)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[P] - current.epsilon_tilde[-K+P+Q] - current.epsilon_tilde[-Q]) 
	* current.interactions_differing_spin(-K, P, P+Q);
one_value -= 8.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[-K-Q] + current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-K] - current.epsilon_tilde[-P]) 
	* current.interactions_differing_spin(-K-Q, -P+Q, Q);
one_value -= 8.000000 * current.interactions_same_spin(K, P, Q)
	* sign(current.epsilon_tilde[-K-Q] + current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-K] - current.epsilon_tilde[-P]) 
	* current.interactions_same_spin(-K-Q, -P+Q, Q);
one_value += 8.000000 * current.interactions_same_spin(-Q, -K+P+Q, P+Q)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[P] - current.epsilon_tilde[-K+P+Q] - current.epsilon_tilde[-Q]) 
	* current.interactions_same_spin(-K, P, P+Q);
nQ_value -= 8.000000 * current.interactions_same_spin(-K, K-P-Q, K-P)
	* sign(current.epsilon_tilde[P] + current.epsilon_tilde[Q] - current.epsilon_tilde[K] - current.epsilon_tilde[-K+P+Q]) 
	* current.interactions_same_spin(P, Q, K-P);
nQ_value += 8.000000 * current.interactions_same_spin(K, P, P-Q)
	* sign(current.epsilon_tilde[Q] + current.epsilon_tilde[K+P-Q] - current.epsilon_tilde[P] - current.epsilon_tilde[K]) 
	* current.interactions_same_spin(Q, K+P-Q, P-Q);
nQ_value -= 8.000000 * current.interactions_same_spin(P, K-P-Q, K-P)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[Q] - current.epsilon_tilde[-P] - current.epsilon_tilde[-K+P+Q]) 
	* current.interactions_same_spin(-K, Q, K-P);
nQ_value += 8.000000 * current.interactions_same_spin(-K, -P, K-Q)
	* sign(current.epsilon_tilde[Q] + current.epsilon_tilde[K+P-Q] - current.epsilon_tilde[K] - current.epsilon_tilde[P]) 
	* current.interactions_same_spin(Q, K+P-Q, K-Q);
nQ_value += 8.000000 * current.interactions_same_spin(P, Q, K-P)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[K-P-Q] - current.epsilon_tilde[-P] - current.epsilon_tilde[-Q]) 
	* current.interactions_same_spin(-K, K-P-Q, K-P);
nQ_value -= 8.000000 * current.interactions_same_spin(-P, -K+P+Q, P+Q)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[Q] - current.epsilon_tilde[-K+P+Q] - current.epsilon_tilde[-P]) 
	* current.interactions_same_spin(-K, Q, P+Q);
nQ_value -= 8.000000 * current.interactions_differing_spin(-K, K-P-Q, K-P)
	* sign(current.epsilon_tilde[P] + current.epsilon_tilde[Q] - current.epsilon_tilde[K] - current.epsilon_tilde[-K+P+Q]) 
	* current.interactions_differing_spin(P, Q, K-P);
nQ_value += 8.000000 * current.interactions_differing_spin(K, P, P-Q)
	* sign(current.epsilon_tilde[Q] + current.epsilon_tilde[K+P-Q] - current.epsilon_tilde[P] - current.epsilon_tilde[K]) 
	* current.interactions_differing_spin(Q, K+P-Q, P-Q);
nQ_value -= 8.000000 * current.interactions_differing_spin(P, K-P-Q, K-P)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[Q] - current.epsilon_tilde[-P] - current.epsilon_tilde[-K+P+Q]) 
	* current.interactions_differing_spin(-K, Q, K-P);
nQ_value += 8.000000 * current.interactions_differing_spin(-K, -P, K-Q)
	* sign(current.epsilon_tilde[Q] + current.epsilon_tilde[K+P-Q] - current.epsilon_tilde[K] - current.epsilon_tilde[P]) 
	* current.interactions_differing_spin(Q, K+P-Q, K-Q);
nQ_value += 8.000000 * current.interactions_differing_spin(P, Q, K-P)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[K-P-Q] - current.epsilon_tilde[-P] - current.epsilon_tilde[-Q]) 
	* current.interactions_differing_spin(-K, K-P-Q, K-P);
nQ_value -= 8.000000 * current.interactions_differing_spin(-P, -K+P+Q, P+Q)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[Q] - current.epsilon_tilde[-K+P+Q] - current.epsilon_tilde[-P]) 
	* current.interactions_differing_spin(-K, Q, P+Q);
nQ_value += 8.000000 * current.interactions_same_spin(P, K-P-Q, K-P)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[Q] - current.epsilon_tilde[-K+P+Q] - current.epsilon_tilde[-P]) 
	* current.interactions_same_spin(-K, Q, P+Q);
nQ_value -= 8.000000 * current.interactions_same_spin(-K, -P, K-Q)
	* sign(current.epsilon_tilde[Q] + current.epsilon_tilde[K+P-Q] - current.epsilon_tilde[P] - current.epsilon_tilde[K]) 
	* current.interactions_same_spin(Q, K+P-Q, P-Q);
nQ_value -= 8.000000 * current.interactions_same_spin(P, Q, K-P)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[K-P-Q] - current.epsilon_tilde[-Q] - current.epsilon_tilde[-P]) 
	* current.interactions_same_spin(-K, K-P-Q, K-Q);
nQ_value -= 8.000000 * current.interactions_same_spin(K, P, P-Q)
	* sign(current.epsilon_tilde[Q] + current.epsilon_tilde[K+P-Q] - current.epsilon_tilde[K] - current.epsilon_tilde[P]) 
	* current.interactions_same_spin(Q, K+P-Q, K-Q);
nQ_value += 8.000000 * current.interactions_same_spin(-P, -K+P+Q, P+Q)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[Q] - current.epsilon_tilde[-P] - current.epsilon_tilde[-K+P+Q]) 
	* current.interactions_same_spin(-K, Q, K-P);
nQ_value += 8.000000 * current.interactions_same_spin(-K, K-P-Q, K-P)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[-Q] - current.epsilon_tilde[K-P-Q] - current.epsilon_tilde[-K]) 
	* current.interactions_same_spin(-P, -Q, K-Q);
nQ_value *= occupation_numbers[Q];
} // r-loop
dHdl.dispersion[K] += (nQ_value + one_value) * occupation_numbers[P];
} // q-loop
} // p-loop

//--------------------------------------------------------------//
#pragma omp parallel for
for (int K_pos=0; K_pos < N; ++K_pos) {
momentum_iterator<L> K(K_pos);
for (momentum_iterator<L> P = momentum_iterator<L>::begin(); P != momentum_iterator<L>::end(); ++P) {
for (momentum_iterator<L> Q = momentum_iterator<L>::begin(); Q != momentum_iterator<L>::end(); ++Q) {
dHdl.interactions_differing_spin(K, P, Q)+= 2.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q)
	* current.epsilon_tilde[K+Q];
dHdl.interactions_differing_spin(K, P, Q)-= 4.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q)
	* current.epsilon_tilde[P];
dHdl.interactions_differing_spin(K, P, Q)+= 2.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q)
	* current.epsilon_tilde[P-Q];
for (momentum_iterator<L> R = momentum_iterator<L>::begin(); R != momentum_iterator<L>::end(); ++R) {
double nR_value{};
double one_value{};
one_value += 4.000000 * current.interactions_differing_spin(-R, K+P+R, K+Q+R)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-P] - current.epsilon_tilde[R] - current.epsilon_tilde[-K-P-R]) 
	* current.interactions_differing_spin(-K, -P, K+R);
one_value -= 4.000000 * current.interactions_differing_spin(K, P, R)
	* sign(current.epsilon_tilde[K+R] + current.epsilon_tilde[P-R] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K+R, P-R, Q-R);
nR_value += 8.000000 * current.interactions_differing_spin(R, P-Q, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value -= 8.000000 * current.interactions_differing_spin(K, R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[Q-R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-R]) 
	* current.interactions_differing_spin(-P, Q-R, Q);
nR_value -= 8.000000 * current.interactions_differing_spin(K, R, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value += 8.000000 * current.interactions_differing_spin(K, Q-R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-Q+R]) 
	* current.interactions_differing_spin(-P, R, Q);
nR_value += 8.000000 * current.interactions_differing_spin(R, P-Q, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value -= 8.000000 * current.interactions_differing_spin(K, R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[Q-R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-R]) 
	* current.interactions_differing_spin(-P, Q-R, Q);
nR_value -= 8.000000 * current.interactions_differing_spin(K, R, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value += 8.000000 * current.interactions_differing_spin(K, Q-R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-Q+R]) 
	* current.interactions_differing_spin(-P, R, Q);
nR_value -= 8.000000 * current.interactions_differing_spin(K, Q-R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[R] - current.epsilon_tilde[-Q+R] - current.epsilon_tilde[-P+Q]) 
	* current.interactions_differing_spin(-P, R, P-Q+R);
nR_value += 8.000000 * current.interactions_differing_spin(-K, -R, K+Q-R)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[Q-R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-R]) 
	* current.interactions_differing_spin(-P, Q-R, Q);
nR_value += 8.000000 * current.interactions_differing_spin(-K, -R, K-R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value -= 4.000000 * current.interactions_differing_spin(-R, K+P+R, K+Q+R)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-P] - current.epsilon_tilde[R] - current.epsilon_tilde[-K-P-R]) 
	* current.interactions_differing_spin(-K, -P, K+R);
nR_value += 4.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[P-Q] - current.epsilon_tilde[P-Q] - current.epsilon_tilde[R]) 
	* current.interactions_differing_spin(R, P-Q, P-Q-R);
nR_value += 4.000000 * current.interactions_differing_spin(K, P, P-R)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+P-R] - current.epsilon_tilde[P-Q] - current.epsilon_tilde[K+Q]) 
	* current.interactions_differing_spin(R, K+P-R, P-Q-R);
nR_value -= 4.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+Q] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[R]) 
	* current.interactions_differing_spin(R, K+Q, K+Q-R);
nR_value -= 8.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[-R] - current.epsilon_tilde[-R] - current.epsilon_tilde[-P]) 
	* current.interactions_differing_spin(-P, -R, P-R);
nR_value += 4.000000 * current.interactions_differing_spin(R, P-Q, P-Q-R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value += 4.000000 * current.interactions_differing_spin(-R, K+P+R, P-Q+R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+P+R] - current.epsilon_tilde[-R]) 
	* current.interactions_differing_spin(K, P, P+R);
nR_value -= 4.000000 * current.interactions_differing_spin(R, K+Q, K+Q-R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value -= 4.000000 * current.interactions_differing_spin(-K, -P, K-R)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+P-R] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(R, K+P-R, K+Q-R);
dHdl.interactions_differing_spin(K, P, Q) += one_value + occupation_numbers[R] * nR_value;
} // s-loop
} // r-loop
} // q-loop
} // p-loop

//--------------------------------------------------------------//
#pragma omp parallel for
for (int K_pos=0; K_pos < N; ++K_pos) {
momentum_iterator<L> K(K_pos);
for (momentum_iterator<L> P = momentum_iterator<L>::begin(); P != momentum_iterator<L>::end(); ++P) {
if (K==P) continue; // Pauli principle
for (momentum_iterator<L> Q = momentum_iterator<L>::begin(); Q != momentum_iterator<L>::end(); ++Q) {
if (K+Q==P-Q) continue; // Pauli principle
dHdl.interactions_same_spin(K, P, Q)+= 2.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q)
	* current.epsilon_tilde[K+Q];
dHdl.interactions_same_spin(K, P, Q)-= 4.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q)
	* current.epsilon_tilde[P];
dHdl.interactions_same_spin(K, P, Q)+= 2.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q)
	* current.epsilon_tilde[P-Q];
for (momentum_iterator<L> R = momentum_iterator<L>::begin(); R != momentum_iterator<L>::end(); ++R) {
double nR_value{};
double one_value{};
one_value += 4.000000 * current.interactions_differing_spin(-R, K+P+R, K+Q+R)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-P] - current.epsilon_tilde[R] - current.epsilon_tilde[-K-P-R]) 
	* current.interactions_differing_spin(-K, -P, K+R);
one_value -= 4.000000 * current.interactions_differing_spin(K, P, R)
	* sign(current.epsilon_tilde[K+R] + current.epsilon_tilde[P-R] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K+R, P-R, Q-R);
nR_value += 8.000000 * current.interactions_differing_spin(R, P-Q, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value -= 8.000000 * current.interactions_differing_spin(K, R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[Q-R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-R]) 
	* current.interactions_differing_spin(-P, Q-R, Q);
nR_value -= 8.000000 * current.interactions_differing_spin(K, R, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value += 8.000000 * current.interactions_differing_spin(K, Q-R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-Q+R]) 
	* current.interactions_differing_spin(-P, R, Q);
nR_value += 8.000000 * current.interactions_differing_spin(R, P-Q, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value -= 8.000000 * current.interactions_differing_spin(K, R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[Q-R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-R]) 
	* current.interactions_differing_spin(-P, Q-R, Q);
nR_value -= 8.000000 * current.interactions_differing_spin(K, R, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value += 8.000000 * current.interactions_differing_spin(K, Q-R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-Q+R]) 
	* current.interactions_differing_spin(-P, R, Q);
nR_value -= 8.000000 * current.interactions_differing_spin(K, Q-R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[R] - current.epsilon_tilde[-Q+R] - current.epsilon_tilde[-P+Q]) 
	* current.interactions_differing_spin(-P, R, P-Q+R);
nR_value += 8.000000 * current.interactions_differing_spin(-K, -R, K+Q-R)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[Q-R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-R]) 
	* current.interactions_differing_spin(-P, Q-R, Q);
nR_value += 8.000000 * current.interactions_differing_spin(-K, -R, K-R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value -= 4.000000 * current.interactions_differing_spin(-R, K+P+R, K+Q+R)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-P] - current.epsilon_tilde[R] - current.epsilon_tilde[-K-P-R]) 
	* current.interactions_differing_spin(-K, -P, K+R);
nR_value += 4.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[P-Q] - current.epsilon_tilde[P-Q] - current.epsilon_tilde[R]) 
	* current.interactions_differing_spin(R, P-Q, P-Q-R);
nR_value += 4.000000 * current.interactions_differing_spin(K, P, P-R)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+P-R] - current.epsilon_tilde[P-Q] - current.epsilon_tilde[K+Q]) 
	* current.interactions_differing_spin(R, K+P-R, P-Q-R);
nR_value -= 4.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+Q] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[R]) 
	* current.interactions_differing_spin(R, K+Q, K+Q-R);
nR_value -= 8.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[-R] - current.epsilon_tilde[-R] - current.epsilon_tilde[-P]) 
	* current.interactions_differing_spin(-P, -R, P-R);
nR_value += 4.000000 * current.interactions_differing_spin(R, P-Q, P-Q-R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value += 4.000000 * current.interactions_differing_spin(-R, K+P+R, P-Q+R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+P+R] - current.epsilon_tilde[-R]) 
	* current.interactions_differing_spin(K, P, P+R);
nR_value -= 4.000000 * current.interactions_differing_spin(R, K+Q, K+Q-R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
nR_value -= 4.000000 * current.interactions_differing_spin(-K, -P, K-R)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+P-R] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(R, K+P-R, K+Q-R);
dHdl.interactions_same_spin(K, P, Q) += one_value + occupation_numbers[R] * nR_value;
} // s-loop
} // r-loop
} // q-loop
} // p-loop
dHdl.interactions_same_spin.symmetrize();
dHdl.interactions_differing_spin.symmetrize();
}
} // namespace NickelCUT::flow