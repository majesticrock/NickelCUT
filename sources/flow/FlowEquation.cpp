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
for (int p_pos=0; p_pos < N; ++p_pos) {
momentum_iterator<L> p(p_pos);
for (momentum_iterator<L> q = momentum_iterator<L>::begin(); q != momentum_iterator<L>::end(); ++q) {
double nr_value{};
double one_value{};
for (momentum_iterator<L> r = momentum_iterator<L>::begin(); r != momentum_iterator<L>::end(); ++r) {
one_value += 8.000000 * current.interactions_same_spin(K, P, Q)
	* sign(current.epsilon_tilde[-K-Q] + current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-P] - current.epsilon_tilde[-K]) 
	* current.interactions_same_spin(-K-Q, -P+Q, K-P+Q);
one_value -= 8.000000 * current.interactions_same_spin(-Q, -K+P+Q, P+Q)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[P] - current.epsilon_tilde[-Q] - current.epsilon_tilde[-K+P+Q]) 
	* current.interactions_same_spin(-K, P, K-Q);
one_value += 8.000000 * current.interactions_differing_spin(-Q, -K+P+Q, P+Q)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[P] - current.epsilon_tilde[-K+P+Q] - current.epsilon_tilde[-Q]) 
	* current.interactions_differing_spin(-K, P, P+Q);
one_value -= 8.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[-K-Q] + current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-K] - current.epsilon_tilde[-P]) 
	* current.interactions_differing_spin(-K-Q, -P+Q, Q);
nr_value += 4.000000 * current.interactions_differing_spin(P, Q, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[Q] - current.epsilon_tilde[K] - current.epsilon_tilde[Q]) 
	* current.interactions_differing_spin(K, Q, Gamma<L>);
nr_value -= 4.000000 * current.interactions_differing_spin(P, Q, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K] - current.epsilon_tilde[P]) 
	* current.interactions_differing_spin(K, P, Gamma<L>);
nr_value -= 4.000000 * current.interactions_same_spin(-P, -Q, P-Q)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[Q] - current.epsilon_tilde[K] - current.epsilon_tilde[Q]) 
	* current.interactions_differing_spin(K, Q, Gamma<L>);
nr_value += 4.000000 * current.interactions_same_spin(-P, -Q, P-Q)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K] - current.epsilon_tilde[P]) 
	* current.interactions_differing_spin(K, P, Gamma<L>);
nr_value -= 8.000000 * current.interactions_differing_spin(-P, -K+P+Q, P+Q)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[Q] - current.epsilon_tilde[-K+P+Q] - current.epsilon_tilde[-P]) 
	* current.interactions_differing_spin(-K, Q, P+Q);
nr_value -= 8.000000 * current.interactions_differing_spin(P, K-P-Q, K-P)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[Q] - current.epsilon_tilde[-P] - current.epsilon_tilde[-K+P+Q]) 
	* current.interactions_differing_spin(-K, Q, K-P);
nr_value -= 8.000000 * current.interactions_differing_spin(-K, K-P-Q, K-P)
	* sign(current.epsilon_tilde[P] + current.epsilon_tilde[Q] - current.epsilon_tilde[K] - current.epsilon_tilde[-K+P+Q]) 
	* current.interactions_differing_spin(P, Q, K-P);
nr_value += 8.000000 * current.interactions_differing_spin(K, P, P-Q)
	* sign(current.epsilon_tilde[Q] + current.epsilon_tilde[K+P-Q] - current.epsilon_tilde[P] - current.epsilon_tilde[K]) 
	* current.interactions_differing_spin(Q, K+P-Q, P-Q);
nr_value += 8.000000 * current.interactions_differing_spin(P, Q, K-P)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[K-P-Q] - current.epsilon_tilde[-P] - current.epsilon_tilde[-Q]) 
	* current.interactions_differing_spin(-K, K-P-Q, K-P);
nr_value += 8.000000 * current.interactions_differing_spin(-K, -P, K-Q)
	* sign(current.epsilon_tilde[Q] + current.epsilon_tilde[K+P-Q] - current.epsilon_tilde[K] - current.epsilon_tilde[P]) 
	* current.interactions_differing_spin(Q, K+P-Q, K-Q);
nr_value += 4.000000 * current.interactions_differing_spin(P, Q, Gamma<L>)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-P] - current.epsilon_tilde[-P] - current.epsilon_tilde[-K]) 
	* current.interactions_same_spin(-K, -P, K-P);
nr_value -= 4.000000 * current.interactions_differing_spin(P, Q, Gamma<L>)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-Q] - current.epsilon_tilde[-Q] - current.epsilon_tilde[-K]) 
	* current.interactions_same_spin(-K, -Q, K-Q);
nr_value -= 8.000000 * current.interactions_same_spin(K, P, P-Q)
	* sign(current.epsilon_tilde[Q] + current.epsilon_tilde[K+P-Q] - current.epsilon_tilde[K] - current.epsilon_tilde[P]) 
	* current.interactions_same_spin(Q, K+P-Q, K-Q);
nr_value += 4.000000 * current.interactions_same_spin(-P, -Q, P-Q)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-Q] - current.epsilon_tilde[-Q] - current.epsilon_tilde[-K]) 
	* current.interactions_same_spin(-K, -Q, K-Q);
nr_value += 8.000000 * current.interactions_same_spin(P, K-P-Q, K-P)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[Q] - current.epsilon_tilde[-K+P+Q] - current.epsilon_tilde[-P]) 
	* current.interactions_same_spin(-K, Q, P+Q);
nr_value += 8.000000 * current.interactions_same_spin(-P, -K+P+Q, P+Q)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[Q] - current.epsilon_tilde[-P] - current.epsilon_tilde[-K+P+Q]) 
	* current.interactions_same_spin(-K, Q, K-P);
nr_value -= 8.000000 * current.interactions_same_spin(-K, -P, K-Q)
	* sign(current.epsilon_tilde[Q] + current.epsilon_tilde[K+P-Q] - current.epsilon_tilde[P] - current.epsilon_tilde[K]) 
	* current.interactions_same_spin(Q, K+P-Q, P-Q);
nr_value -= 8.000000 * current.interactions_same_spin(P, Q, K-P)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[K-P-Q] - current.epsilon_tilde[-Q] - current.epsilon_tilde[-P]) 
	* current.interactions_same_spin(-K, K-P-Q, K-Q);
nr_value -= 4.000000 * current.interactions_same_spin(-P, -Q, P-Q)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-P] - current.epsilon_tilde[-P] - current.epsilon_tilde[-K]) 
	* current.interactions_same_spin(-K, -P, K-P);
nr_value += 8.000000 * current.interactions_same_spin(-K, K-P-Q, K-P)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[-Q] - current.epsilon_tilde[K-P-Q] - current.epsilon_tilde[-K]) 
	* current.interactions_same_spin(-P, -Q, K-Q);
nr_value *= occupation_numbers[r];
} // r-loop
dHdl.dispersion[p] += (nr_value + one_value) * occupation_numbers[q];
} // q-loop
} // p-loop

//--------------------------------------------------------------//
#pragma omp parallel for
for (int p_pos=0; p_pos < N; ++p_pos) {
momentum_iterator<L> p(p_pos);
for (momentum_iterator<L> q = momentum_iterator<L>::begin(); q != momentum_iterator<L>::end(); ++q) {
for (momentum_iterator<L> r = momentum_iterator<L>::begin(); r != momentum_iterator<L>::end(); ++r) {
dHdl.interactions_differing_spin(p, q, r)+= 4.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q)
	* current.epsilon_tilde[P-Q];
dHdl.interactions_differing_spin(p, q, r)-= 4.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q)
	* current.epsilon_tilde[K];
for (momentum_iterator<L> s = momentum_iterator<L>::begin(); s != momentum_iterator<L>::end(); ++s) {
double ns_value{};
double one_value{};
one_value += 4.000000 * current.interactions_differing_spin(-R, K+P+R, K+Q+R)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-P] - current.epsilon_tilde[R] - current.epsilon_tilde[-K-P-R]) 
	* current.interactions_differing_spin(-K, -P, K+R);
one_value -= 4.000000 * current.interactions_differing_spin(K, P, R)
	* sign(current.epsilon_tilde[K+R] + current.epsilon_tilde[P-R] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K+R, P-R, Q-R);
ns_value -= 4.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+Q] - current.epsilon_tilde[R] - current.epsilon_tilde[K+Q]) 
	* current.interactions_differing_spin(R, K+Q, Gamma<L>);
ns_value -= 4.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[P-Q] - current.epsilon_tilde[R] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(R, P-Q, Gamma<L>);
ns_value += 8.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[P] + current.epsilon_tilde[R] - current.epsilon_tilde[P] - current.epsilon_tilde[R]) 
	* current.interactions_same_spin(P, R, Gamma<L>);
ns_value += 4.000000 * current.interactions_same_spin(R, P-Q, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
ns_value += 4.000000 * current.interactions_differing_spin(R, K+Q, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
ns_value += 8.000000 * current.interactions_differing_spin(K, Q-R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-Q+R]) 
	* current.interactions_same_spin(-P, R, Q);
ns_value -= 8.000000 * current.interactions_differing_spin(K, R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[Q-R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-R]) 
	* current.interactions_same_spin(-P, Q-R, Q);
ns_value -= 8.000000 * current.interactions_differing_spin(K, R, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
ns_value += 4.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+Q] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[R]) 
	* current.interactions_same_spin(R, K+Q, K+Q-R);
ns_value -= 4.000000 * current.interactions_same_spin(R, P-Q, P-Q-R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
ns_value -= 4.000000 * current.interactions_differing_spin(-R, K+P+R, K+Q+R)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-P] - current.epsilon_tilde[R] - current.epsilon_tilde[-K-P-R]) 
	* current.interactions_differing_spin(-K, -P, K+R);
ns_value -= 4.000000 * current.interactions_same_spin(R, K+Q, K+Q-R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
ns_value -= 8.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[-R] - current.epsilon_tilde[-R] - current.epsilon_tilde[-P]) 
	* current.interactions_same_spin(-P, -R, P-R);
ns_value -= 8.000000 * current.interactions_differing_spin(-K, K-P+Q-R, K-R)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[-R] - current.epsilon_tilde[K-P+Q-R] - current.epsilon_tilde[-K-Q]) 
	* current.interactions_differing_spin(-P, -R, K+Q-R);
ns_value -= 8.000000 * current.interactions_differing_spin(K, Q-R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[R] - current.epsilon_tilde[-Q+R] - current.epsilon_tilde[-P+Q]) 
	* current.interactions_same_spin(-P, R, P-Q+R);
ns_value += 8.000000 * current.interactions_differing_spin(-K, -R, P-Q-R)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[-K+P-Q-R] - current.epsilon_tilde[-R] - current.epsilon_tilde[-K-Q]) 
	* current.interactions_differing_spin(-P, -K+P-Q-R, P-R);
ns_value -= 4.000000 * current.interactions_differing_spin(-R, K+P+R, P-Q+R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+P+R] - current.epsilon_tilde[-R]) 
	* current.interactions_differing_spin(K, P, P+R);
ns_value += 4.000000 * current.interactions_differing_spin(K, P, P-R)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+P-R] - current.epsilon_tilde[P-Q] - current.epsilon_tilde[K+Q]) 
	* current.interactions_differing_spin(R, K+P-R, P-Q-R);
ns_value += 4.000000 * current.interactions_differing_spin(-K, -P, K-R)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+P-R] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(R, K+P-R, K+Q-R);
ns_value += 4.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[P-Q] - current.epsilon_tilde[P-Q] - current.epsilon_tilde[R]) 
	* current.interactions_same_spin(R, P-Q, P-Q-R);
ns_value -= 8.000000 * current.interactions_same_spin(-K, -Q+R, K+R)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-Q+R]) 
	* current.interactions_differing_spin(-P, R, Q);
ns_value += 8.000000 * current.interactions_same_spin(-K, -R, K+Q-R)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[Q-R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-R]) 
	* current.interactions_differing_spin(-P, Q-R, Q);
ns_value += 8.000000 * current.interactions_differing_spin(K, R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[Q-R] - current.epsilon_tilde[-R] - current.epsilon_tilde[-P+Q]) 
	* current.interactions_same_spin(-P, Q-R, P-R);
ns_value += 8.000000 * current.interactions_same_spin(-K, -R, K-R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q);
dHdl.interactions_differing_spin(p, q, r) += one_value + occupation_numbers[s] * ns_value;
} // s-loop
} // r-loop
} // q-loop
} // p-loop

//--------------------------------------------------------------//
#pragma omp parallel for
for (int p_pos=0; p_pos < N; ++p_pos) {
momentum_iterator<L> p(p_pos);
for (momentum_iterator<L> q = momentum_iterator<L>::begin(); q != momentum_iterator<L>::end(); ++q) {
if (p==q) continue; // Pauli principle
for (momentum_iterator<L> r = momentum_iterator<L>::begin(); r != momentum_iterator<L>::end(); ++r) {
if (p+r==q-r) continue; // Pauli principle
dHdl.interactions_same_spin(p, q, r)+= 4.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(K, P, Q)
	* current.epsilon_tilde[P-Q];
dHdl.interactions_same_spin(p, q, r)-= 4.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(K, P, Q)
	* current.epsilon_tilde[K];
for (momentum_iterator<L> s = momentum_iterator<L>::begin(); s != momentum_iterator<L>::end(); ++s) {
double ns_value{};
double one_value{};
one_value += 4.000000 * current.interactions_same_spin(-R, K+P+R, K+Q+R)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-P] - current.epsilon_tilde[R] - current.epsilon_tilde[-K-P-R]) 
	* current.interactions_same_spin(-K, -P, K+R);
one_value -= 4.000000 * current.interactions_same_spin(K, P, R)
	* sign(current.epsilon_tilde[K+R] + current.epsilon_tilde[P-R] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(K+R, P-R, Q-R);
ns_value -= 4.000000 * current.interactions_same_spin(K, P, Q)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+Q] - current.epsilon_tilde[R] - current.epsilon_tilde[K+Q]) 
	* current.interactions_differing_spin(R, K+Q, Gamma<L>);
ns_value -= 4.000000 * current.interactions_same_spin(K, P, Q)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[P-Q] - current.epsilon_tilde[R] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(R, P-Q, Gamma<L>);
ns_value += 8.000000 * current.interactions_same_spin(K, P, Q)
	* sign(current.epsilon_tilde[P] + current.epsilon_tilde[R] - current.epsilon_tilde[P] - current.epsilon_tilde[R]) 
	* current.interactions_differing_spin(P, R, Gamma<L>);
ns_value += 4.000000 * current.interactions_differing_spin(R, P-Q, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(K, P, Q);
ns_value += 4.000000 * current.interactions_differing_spin(R, K+Q, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(K, P, Q);
ns_value += 8.000000 * current.interactions_differing_spin(K, Q-R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-Q+R]) 
	* current.interactions_differing_spin(-P, R, Q);
ns_value -= 8.000000 * current.interactions_differing_spin(K, R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[Q-R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-R]) 
	* current.interactions_differing_spin(-P, Q-R, Q);
ns_value -= 8.000000 * current.interactions_differing_spin(K, R, Gamma<L>)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(K, P, Q);
ns_value += 4.000000 * current.interactions_same_spin(K, P, Q)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+Q] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[R]) 
	* current.interactions_same_spin(R, K+Q, K+Q-R);
ns_value -= 4.000000 * current.interactions_same_spin(R, P-Q, P-Q-R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(K, P, Q);
ns_value -= 4.000000 * current.interactions_same_spin(-R, K+P+R, K+Q+R)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-P] - current.epsilon_tilde[R] - current.epsilon_tilde[-K-P-R]) 
	* current.interactions_same_spin(-K, -P, K+R);
ns_value -= 4.000000 * current.interactions_same_spin(R, K+Q, K+Q-R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(K, P, Q);
ns_value -= 8.000000 * current.interactions_same_spin(K, P, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[-R] - current.epsilon_tilde[-R] - current.epsilon_tilde[-P]) 
	* current.interactions_same_spin(-P, -R, P-R);
ns_value -= 8.000000 * current.interactions_same_spin(-K, K-P+Q-R, K-R)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[-R] - current.epsilon_tilde[K-P+Q-R] - current.epsilon_tilde[-K-Q]) 
	* current.interactions_same_spin(-P, -R, K+Q-R);
ns_value -= 8.000000 * current.interactions_same_spin(K, Q-R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[R] - current.epsilon_tilde[-Q+R] - current.epsilon_tilde[-P+Q]) 
	* current.interactions_same_spin(-P, R, P-Q+R);
ns_value += 8.000000 * current.interactions_same_spin(-K, -R, P-Q-R)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[-K+P-Q-R] - current.epsilon_tilde[-R] - current.epsilon_tilde[-K-Q]) 
	* current.interactions_same_spin(-P, -K+P-Q-R, P-R);
ns_value -= 4.000000 * current.interactions_same_spin(-R, K+P+R, P-Q+R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+P+R] - current.epsilon_tilde[-R]) 
	* current.interactions_same_spin(K, P, P+R);
ns_value += 4.000000 * current.interactions_same_spin(K, P, P-R)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+P-R] - current.epsilon_tilde[P-Q] - current.epsilon_tilde[K+Q]) 
	* current.interactions_same_spin(R, K+P-R, P-Q-R);
ns_value += 4.000000 * current.interactions_same_spin(-K, -P, K-R)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[K+P-R] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(R, K+P-R, K+Q-R);
ns_value += 4.000000 * current.interactions_same_spin(K, P, Q)
	* sign(current.epsilon_tilde[R] + current.epsilon_tilde[P-Q] - current.epsilon_tilde[P-Q] - current.epsilon_tilde[R]) 
	* current.interactions_same_spin(R, P-Q, P-Q-R);
ns_value -= 8.000000 * current.interactions_same_spin(-K, -Q+R, K+R)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-Q+R]) 
	* current.interactions_same_spin(-P, R, Q);
ns_value += 8.000000 * current.interactions_same_spin(-K, -R, K+Q-R)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[Q-R] - current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-R]) 
	* current.interactions_same_spin(-P, Q-R, Q);
ns_value += 8.000000 * current.interactions_same_spin(K, R, Q)
	* sign(current.epsilon_tilde[-P] + current.epsilon_tilde[Q-R] - current.epsilon_tilde[-R] - current.epsilon_tilde[-P+Q]) 
	* current.interactions_same_spin(-P, Q-R, P-R);
ns_value += 8.000000 * current.interactions_same_spin(-K, -R, K-R)
	* sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(K, P, Q);
dHdl.interactions_same_spin(p, q, r) += one_value + occupation_numbers[s] * ns_value;
} // s-loop
} // r-loop
} // q-loop
} // p-loop
dHdl.interactions_same_spin.symmetrize();
dHdl.interactions_differing_spin.symmetrize();
}
} // namespace NickelCUT::flow