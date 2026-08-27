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
one_value -= 8.000000 * current.interactions_differing_spin(K, P, Q)
	* sign(current.epsilon_tilde[-K-Q] + current.epsilon_tilde[-P+Q] - current.epsilon_tilde[-K] - current.epsilon_tilde[-P]) 
	* current.interactions_differing_spin(-K-Q, -P+Q, Q);
one_value += 8.000000 * current.interactions_differing_spin(-Q, -K+P+Q, P+Q)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[P] - current.epsilon_tilde[-K+P+Q] - current.epsilon_tilde[-Q]) 
	* current.interactions_differing_spin(-K, P, P+Q);
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
dHdl.interactions_differing_spin(K, P, Q)-= 4.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q)
	* current.epsilon_tilde[K];
dHdl.interactions_differing_spin(K, P, Q)+= 4.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K, P, Q)
	* current.epsilon_tilde[P-Q];
for (momentum_iterator<L> R = momentum_iterator<L>::begin(); R != momentum_iterator<L>::end(); ++R) {
double nR_value{};
double one_value{};
one_value -= 4.000000 * current.interactions_differing_spin(K, P, R)
	* sign(current.epsilon_tilde[K+R] + current.epsilon_tilde[P-R] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_differing_spin(K+R, P-R, Q-R);
one_value += 4.000000 * current.interactions_differing_spin(-R, K+P+R, K+Q+R)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-P] - current.epsilon_tilde[R] - current.epsilon_tilde[-K-P-R]) 
	* current.interactions_differing_spin(-K, -P, K+R);
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
dHdl.interactions_same_spin(K, P, Q)-= 4.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(K, P, Q)
	* current.epsilon_tilde[K];
dHdl.interactions_same_spin(K, P, Q)+= 4.000000 * sign(current.epsilon_tilde[K] + current.epsilon_tilde[P] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(K, P, Q)
	* current.epsilon_tilde[P-Q];
for (momentum_iterator<L> R = momentum_iterator<L>::begin(); R != momentum_iterator<L>::end(); ++R) {
double nR_value{};
double one_value{};
one_value -= 4.000000 * current.interactions_same_spin(K, P, R)
	* sign(current.epsilon_tilde[K+R] + current.epsilon_tilde[P-R] - current.epsilon_tilde[K+Q] - current.epsilon_tilde[P-Q]) 
	* current.interactions_same_spin(K+R, P-R, Q-R);
one_value += 4.000000 * current.interactions_same_spin(-R, K+P+R, K+Q+R)
	* sign(current.epsilon_tilde[-K] + current.epsilon_tilde[-P] - current.epsilon_tilde[R] - current.epsilon_tilde[-K-P-R]) 
	* current.interactions_same_spin(-K, -P, K+R);
dHdl.interactions_same_spin(K, P, Q) += one_value + occupation_numbers[R] * nR_value;
} // s-loop
} // r-loop
} // q-loop
} // p-loop
dHdl.interactions_same_spin.symmetrize();
dHdl.interactions_differing_spin.symmetrize();

static int i=0; ++i;

for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
    if (!NickelCUT::float_equal(dHdl.dispersion[p], dHdl.dispersion[-p])) {
        std::cerr << i << ": Dispersion is not inversion symmetric: " << dHdl.dispersion[p] << "  " << dHdl.dispersion[-p] << std::endl;
        break;
    }
}

bool error_found = false;
for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end() && !error_found; ++p) {
    for (momentum_iterator<L> q = momentum_iterator<L>::begin(); q != momentum_iterator<L>::end() && !error_found; ++q) {
        for (momentum_iterator<L> r = momentum_iterator<L>::begin(); r != momentum_iterator<L>::end() && !error_found; ++r) {
            if(!NickelCUT::float_equal(dHdl.interactions_differing_spin(p, q, r), dHdl.interactions_differing_spin(-p, -q, -r))) {
                std::cerr << i << ": Interaction is not inversion symmetric " 
                    << dHdl.interactions_differing_spin(p, q, r) 
                    << "  " << dHdl.interactions_differing_spin(-p, -q, -r)
                    << std::endl;
                error_found = true;
            }
        }
    }
}
error_found = false;
for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end() && !error_found; ++p) {
    for (momentum_iterator<L> q = momentum_iterator<L>::begin(); q != momentum_iterator<L>::end() && !error_found; ++q) {
        for (momentum_iterator<L> r = momentum_iterator<L>::begin(); r != momentum_iterator<L>::end() && !error_found; ++r) {
            if(!NickelCUT::float_equal(dHdl.interactions_differing_spin(p, q, r), dHdl.interactions_differing_spin(p+r, q-r, -r))) {
                std::cerr << i << ": Interaction is not Hermitian " 
                    << dHdl.interactions_differing_spin(p, q, r) 
                    << "  " << dHdl.interactions_differing_spin(p+r, q-r, -r)
                    << std::endl;
                error_found = true;
            }
        }
    }
}
}
} // namespace NickelCUT::flow