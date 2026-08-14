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
one_value -= 4.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p-r]) 
	* current.interactions_same_spin(-p, -q, p-q+r);
one_value += 4.000000 * current.interactions_differing_spin(-r, -p+q+r, q+r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[q] - current.epsilon_tilde[-p+q+r] - current.epsilon_tilde[-r]) 
	* current.interactions_differing_spin(-p, q, q+r);
one_value -= 4.000000 * current.interactions_differing_spin(-r, -p+q+r, q+r)
	* sign(current.epsilon_tilde[-r] + current.epsilon_tilde[-p+q+r] - current.epsilon_tilde[q] - current.epsilon_tilde[-p]) 
	* current.interactions_differing_spin(-r, -p+q+r, q+r);
one_value -= 4.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, -q+r, r);
one_value += 4.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r);
one_value += 4.000000 * current.interactions_same_spin(-r, -p+q+r, q+r)
	* sign(current.epsilon_tilde[r] + current.epsilon_tilde[p-q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(r, p-q-r, p-r);
one_value += 4.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-q] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p-r, -q+r, p-q+r);
one_value -= 4.000000 * current.interactions_same_spin(-r, -p+q+r, q+r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[q] - current.epsilon_tilde[-r] - current.epsilon_tilde[-p+q+r]) 
	* current.interactions_same_spin(-p, q, p-r);
nr_value += 4.000000 * current.interactions_differing_spin(-p, -r, p-q)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[p-q+r] - current.epsilon_tilde[p] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(q, p-q+r, p-q);
nr_value += 4.000000 * current.interactions_same_spin(p, q, q-r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[-r] - current.epsilon_tilde[-p-q+r]) 
	* current.interactions_same_spin(-p, -q, p-r);
nr_value += 4.000000 * current.interactions_differing_spin(p, q, q-r)
	* sign(current.epsilon_tilde[r] + current.epsilon_tilde[p+q-r] - current.epsilon_tilde[q] - current.epsilon_tilde[p]) 
	* current.interactions_differing_spin(r, p+q-r, q-r);
nr_value -= 4.000000 * current.interactions_differing_spin(q, r, p-q)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[-p+q+r]) 
	* current.interactions_differing_spin(q, r, p-q);
nr_value += 4.000000 * current.interactions_differing_spin(q, r, p-q)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q-r] - current.epsilon_tilde[-q] - current.epsilon_tilde[-r]) 
	* current.interactions_differing_spin(-p, p-q-r, p-q);
nr_value -= 4.000000 * current.interactions_differing_spin(p, q, q-r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+q-r] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(p, q, q-r);
nr_value += 4.000000 * current.interactions_differing_spin(q, p-q-r, p-q)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[p-q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[-r]) 
	* current.interactions_differing_spin(q, p-q-r, p-q);
nr_value += 4.000000 * current.interactions_differing_spin(-q, -p+q+r, q+r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-p+q+r] - current.epsilon_tilde[r] - current.epsilon_tilde[-p]) 
	* current.interactions_differing_spin(-q, -p+q+r, q+r);
nr_value += 4.000000 * current.interactions_differing_spin(-p, p-q-r, p-q)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q-r] - current.epsilon_tilde[-q] - current.epsilon_tilde[-r]) 
	* current.interactions_differing_spin(-p, p-q-r, p-q);
nr_value -= 4.000000 * current.interactions_differing_spin(q, p-q-r, p-q)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[r] - current.epsilon_tilde[-q] - current.epsilon_tilde[-p+q+r]) 
	* current.interactions_differing_spin(-p, r, p-q);
nr_value -= 4.000000 * current.interactions_differing_spin(-q, -p+q+r, q+r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[r] - current.epsilon_tilde[-p+q+r] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p, r, q+r);
nr_value -= 4.000000 * current.interactions_differing_spin(-p, p-q-r, p-q)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[-p+q+r]) 
	* current.interactions_differing_spin(q, r, p-q);
nr_value -= 4.000000 * current.interactions_differing_spin(-p, -r, p-q)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r] - current.epsilon_tilde[-q] - current.epsilon_tilde[-p+q-r]) 
	* current.interactions_differing_spin(-p, -r, p-q);
nr_value -= 4.000000 * current.interactions_same_spin(q, r, p-q)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q-r] - current.epsilon_tilde[-r] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p, p-q-r, p-r);
nr_value += 4.000000 * current.interactions_same_spin(-p, p-q-r, p-q)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-r] - current.epsilon_tilde[p-q-r] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-q, -r, p-r);
nr_value += 4.000000 * current.interactions_same_spin(q, p-q-r, p-q)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[r] - current.epsilon_tilde[-p+q+r] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p, r, q+r);
nr_value += 4.000000 * current.interactions_same_spin(-q, -p+q+r, q+r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[r] - current.epsilon_tilde[-q] - current.epsilon_tilde[-p+q+r]) 
	* current.interactions_same_spin(-p, r, p-q);
nr_value -= 4.000000 * current.interactions_same_spin(-p, p-q-r, p-q)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q-r] - current.epsilon_tilde[-r] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p, p-q-r, p-r);
nr_value -= 4.000000 * current.interactions_same_spin(q, p-q-r, p-q)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-p+q+r] - current.epsilon_tilde[r] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-q, -p+q+r, q+r);
nr_value -= 4.000000 * current.interactions_same_spin(-q, -p+q+r, q+r)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[p-q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[-r]) 
	* current.interactions_same_spin(q, p-q-r, p-q);
nr_value += 4.000000 * current.interactions_same_spin(-p, -q, p-r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+q-r] - current.epsilon_tilde[r]) 
	* current.interactions_same_spin(p, q, q-r);
nr_value -= 4.000000 * current.interactions_same_spin(p, q, q-r)
	* sign(current.epsilon_tilde[r] + current.epsilon_tilde[p+q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(r, p+q-r, p-r);
nr_value += 4.000000 * current.interactions_same_spin(q, r, p-q)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-r] - current.epsilon_tilde[p-q-r] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-q, -r, p-r);
nr_value -= 4.000000 * current.interactions_same_spin(-p, -q, p-r)
	* sign(current.epsilon_tilde[r] + current.epsilon_tilde[p+q-r] - current.epsilon_tilde[q] - current.epsilon_tilde[p]) 
	* current.interactions_same_spin(r, p+q-r, q-r);
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
dHdl.interactions_differing_spin(p, q, r)+= sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* current.epsilon_tilde[q-r];
dHdl.interactions_differing_spin(p, q, r)-= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* current.epsilon_tilde[q];
dHdl.interactions_differing_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p+r] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_differing_spin(p+r, q-r, -r)
	* current.epsilon_tilde[q];
dHdl.interactions_differing_spin(p, q, r)-= sign(current.epsilon_tilde[p+r] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_differing_spin(p+r, q-r, -r)
	* current.epsilon_tilde[q-r];
dHdl.interactions_differing_spin(p, q, r)-= sign(current.epsilon_tilde[p+r] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_differing_spin(p+r, q-r, -r)
	* current.epsilon_tilde[p+r];
dHdl.interactions_differing_spin(p, q, r)+= sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* current.epsilon_tilde[p+r];
dHdl.interactions_differing_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p+r] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_differing_spin(p+r, q-r, -r)
	* current.epsilon_tilde[p]
	* occupation_numbers[p];
dHdl.interactions_differing_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* current.epsilon_tilde[q]
	* occupation_numbers[q];
dHdl.interactions_differing_spin(p, q, r)-= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* current.epsilon_tilde[p]
	* occupation_numbers[p];
dHdl.interactions_differing_spin(p, q, r)-= 2.000000 * sign(current.epsilon_tilde[p+r] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_differing_spin(p+r, q-r, -r)
	* current.epsilon_tilde[q]
	* occupation_numbers[q];
for (momentum_iterator<L> s = momentum_iterator<L>::begin(); s != momentum_iterator<L>::end(); ++s) {
double ns_value{};
double one_value{};
one_value += 2.000000 * current.interactions_differing_spin(p, q, s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p-s] - current.epsilon_tilde[-q+s]) 
	* current.interactions_differing_spin(-p-r, -q+r, r-s);
one_value -= 2.000000 * current.interactions_differing_spin(-s, p+q+s, p+r+s)
	* sign(current.epsilon_tilde[-s] + current.epsilon_tilde[p+q+s] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_differing_spin(-s, p+q+s, p+s);
one_value += 2.000000 * current.interactions_differing_spin(-s, p+q+s, p+r+s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[s] - current.epsilon_tilde[-p-q-s]) 
	* current.interactions_differing_spin(-p, -q, p+s);
one_value -= 2.000000 * current.interactions_differing_spin(p, q, s)
	* sign(current.epsilon_tilde[p+s] + current.epsilon_tilde[q-s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p+s, q-s, r-s);
ns_value -= 2.000000 * current.interactions_differing_spin(p, q, q-s)
	* sign(current.epsilon_tilde[p+r] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p+q-s] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p+r, q-r, q-r-s);
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r);
ns_value += 2.000000 * current.interactions_differing_spin(-p, -q, p-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+q-s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(s, p+q-s, p+r-s);
ns_value -= 2.000000 * current.interactions_differing_spin(-p, -q, p-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p-q+s]) 
	* current.interactions_differing_spin(-p-r, -q+r, p+r-s);
ns_value += 2.000000 * current.interactions_differing_spin(p, q, q-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+q-s] - current.epsilon_tilde[q-r] - current.epsilon_tilde[p+r]) 
	* current.interactions_differing_spin(s, p+q-s, q-r-s);
ns_value += 4.000000 * current.interactions_differing_spin(-q, -s, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_same_spin(-p, -r-s, p-s);
ns_value += 4.000000 * current.interactions_differing_spin(-q, -s, p+r-s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_differing_spin(-p, p-q+r-s, p-s);
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, q-r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, r+s, r);
ns_value -= 2.000000 * current.interactions_same_spin(s, q-r, GammaPoint<L>)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, -q+r, r);
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[r+s] - current.epsilon_tilde[p]) 
	* current.interactions_differing_spin(s, p+r, r);
ns_value -= 2.000000 * current.interactions_differing_spin(s, p+r, GammaPoint<L>)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, -q+r, r);
ns_value -= 4.000000 * current.interactions_same_spin(q, s, GammaPoint<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r);
ns_value -= 4.000000 * current.interactions_differing_spin(p, r-s, r)
	* sign(current.epsilon_tilde[q-r] + current.epsilon_tilde[r-s] - current.epsilon_tilde[q] - current.epsilon_tilde[-s]) 
	* current.interactions_same_spin(q-r, r-s, r);
ns_value += 4.000000 * current.interactions_differing_spin(p, r-s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-r+s]) 
	* current.interactions_same_spin(-q, s, r);
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, r+s, r);
ns_value += 4.000000 * current.interactions_same_spin(q, s, GammaPoint<L>)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, -q+r, r);
ns_value += 2.000000 * current.interactions_differing_spin(s, p+r, GammaPoint<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r);
ns_value += 2.000000 * current.interactions_same_spin(s, q-r, GammaPoint<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r);
ns_value -= 4.000000 * current.interactions_differing_spin(-q, -s, p+r-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[-p+q-r+s]) 
	* current.interactions_differing_spin(s, q-r, p-s);
ns_value += 4.000000 * current.interactions_differing_spin(-p, p-q+r-s, p-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, p-q+r-s, p+r-s);
ns_value += 4.000000 * current.interactions_same_spin(-p, -r+s, p+s)
	* sign(current.epsilon_tilde[q-r] + current.epsilon_tilde[r-s] - current.epsilon_tilde[q] - current.epsilon_tilde[-s]) 
	* current.interactions_differing_spin(q-r, r-s, r);
ns_value += 4.000000 * current.interactions_differing_spin(p, r-s, r)
	* sign(current.epsilon_tilde[-q+r] + current.epsilon_tilde[-r+s] - current.epsilon_tilde[s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q+r, -r+s, q-r+s);
ns_value -= 2.000000 * current.interactions_differing_spin(-s, p+q+s, p+r+s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[s] - current.epsilon_tilde[-p-q-s]) 
	* current.interactions_differing_spin(-p, -q, p+s);
ns_value -= 2.000000 * current.interactions_same_spin(s, q-r, q-r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r);
ns_value -= 2.000000 * current.interactions_same_spin(s, p+r, p+r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r);
ns_value -= 2.000000 * current.interactions_differing_spin(-s, p+q+s, q-r+s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+q+s] - current.epsilon_tilde[-s]) 
	* current.interactions_differing_spin(p, q, q+s);
ns_value -= 4.000000 * current.interactions_same_spin(-p, -r+s, p+s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-r+s]) 
	* current.interactions_differing_spin(-q, s, r);
ns_value -= 4.000000 * current.interactions_differing_spin(-p, p-q+r-s, p-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-s] - current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_differing_spin(-q, -s, p+r-s);
ns_value -= 4.000000 * current.interactions_differing_spin(p, r-s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[-r+s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_same_spin(-q, s, q-r+s);
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, q-r-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[r+s] - current.epsilon_tilde[p]) 
	* current.interactions_differing_spin(s, p+r, r);
ns_value += 2.000000 * current.interactions_differing_spin(-s, p+q+s, q-r+s)
	* sign(current.epsilon_tilde[-s] + current.epsilon_tilde[p+q+s] - current.epsilon_tilde[q] - current.epsilon_tilde[p]) 
	* current.interactions_differing_spin(-s, p+q+s, q+s);
ns_value += 2.000000 * current.interactions_same_spin(s, p+r, p+r-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, -q+r, r);
ns_value += 2.000000 * current.interactions_same_spin(s, q-r, q-r-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, -q+r, r);
ns_value += 2.000000 * current.interactions_differing_spin(-s, p+q+s, p+r+s)
	* sign(current.epsilon_tilde[-s] + current.epsilon_tilde[p+q+s] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_differing_spin(-s, p+q+s, p+s);
ns_value -= 4.000000 * current.interactions_differing_spin(-q, -s, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[p] - current.epsilon_tilde[r+s]) 
	* current.interactions_same_spin(s, p+r, p-s);
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, -q+r, r);
ns_value -= 4.000000 * current.interactions_same_spin(q, s, GammaPoint<L>)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, -q+r, r)
	* occupation_numbers[q];
ns_value -= 4.000000 * current.interactions_differing_spin(p, s, GammaPoint<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* occupation_numbers[p];
ns_value -= 4.000000 * current.interactions_differing_spin(p, s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r-s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-s]) 
	* current.interactions_same_spin(-q, r-s, r)
	* occupation_numbers[r-s];
ns_value += 4.000000 * current.interactions_differing_spin(p, s, GammaPoint<L>)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, -q+r, r)
	* occupation_numbers[p];
ns_value += 4.000000 * current.interactions_differing_spin(p, s, r)
	* sign(current.epsilon_tilde[-s] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[r-s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-s, -q+r, r)
	* occupation_numbers[r-s];
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[r+s] - current.epsilon_tilde[p]) 
	* current.interactions_differing_spin(s, p+r, r)
	* occupation_numbers[r+s];
ns_value += 4.000000 * current.interactions_same_spin(q, s, GammaPoint<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* occupation_numbers[q];
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, r+s, r)
	* occupation_numbers[r+s];
ns_value -= 4.000000 * current.interactions_same_spin(-p, -s, p-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, -q+r, r)
	* occupation_numbers[p];
ns_value -= 4.000000 * current.interactions_differing_spin(-p, p-q+r-s, p-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, p-q+r-s, p+r-s)
	* occupation_numbers[-p+q-r+s];
ns_value += 4.000000 * current.interactions_same_spin(-p, -s, p+r-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r-s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-s]) 
	* current.interactions_differing_spin(-q, r-s, r)
	* occupation_numbers[r-s];
ns_value -= 4.000000 * current.interactions_differing_spin(p, s, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[q] - current.epsilon_tilde[-r+s]) 
	* current.interactions_same_spin(s, q-r, q-s)
	* occupation_numbers[r-s];
ns_value -= 4.000000 * current.interactions_same_spin(-p, -s, p+r-s)
	* sign(current.epsilon_tilde[-s] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[r-s] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-s, -q+r, r)
	* occupation_numbers[r-s];
ns_value += 4.000000 * current.interactions_differing_spin(-q, -s, p+r-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[-p+q-r+s]) 
	* current.interactions_differing_spin(s, q-r, p-s)
	* occupation_numbers[p-q+r-s];
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, q-r-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[r+s] - current.epsilon_tilde[p]) 
	* current.interactions_differing_spin(s, p+r, r)
	* occupation_numbers[r+s];
ns_value += 4.000000 * current.interactions_same_spin(-p, -s, p-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* occupation_numbers[p];
ns_value += 4.000000 * current.interactions_differing_spin(-p, p-q+r-s, p-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-s] - current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_differing_spin(-q, -s, p+r-s)
	* occupation_numbers[-p+q-r+s];
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, -q+r, r)
	* occupation_numbers[q];
ns_value -= 4.000000 * current.interactions_differing_spin(-q, -s, p+r-s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_differing_spin(-p, p-q+r-s, p-s)
	* occupation_numbers[p-q+r-s];
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, q-r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, r+s, r)
	* occupation_numbers[r+s];
ns_value += 4.000000 * current.interactions_differing_spin(p, s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_same_spin(-q, r-s, q-s)
	* occupation_numbers[r-s];
ns_value -= 4.000000 * current.interactions_differing_spin(-q, -s, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_same_spin(-p, -r-s, p-s)
	* occupation_numbers[r+s];
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* occupation_numbers[q];
ns_value += 4.000000 * current.interactions_differing_spin(-q, -s, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[p] - current.epsilon_tilde[r+s]) 
	* current.interactions_same_spin(s, p+r, p-s)
	* occupation_numbers[r+s];
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
dHdl.interactions_same_spin(p, q, r)+= sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* current.epsilon_tilde[q-r];
dHdl.interactions_same_spin(p, q, r)-= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* current.epsilon_tilde[q];
dHdl.interactions_same_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p+r] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(p+r, q-r, -r)
	* current.epsilon_tilde[q];
dHdl.interactions_same_spin(p, q, r)-= sign(current.epsilon_tilde[p+r] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(p+r, q-r, -r)
	* current.epsilon_tilde[q-r];
dHdl.interactions_same_spin(p, q, r)-= sign(current.epsilon_tilde[p+r] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(p+r, q-r, -r)
	* current.epsilon_tilde[p+r];
dHdl.interactions_same_spin(p, q, r)+= sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* current.epsilon_tilde[p+r];
dHdl.interactions_same_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p+r] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(p+r, q-r, -r)
	* current.epsilon_tilde[p]
	* occupation_numbers[p];
dHdl.interactions_same_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* current.epsilon_tilde[q]
	* occupation_numbers[q];
dHdl.interactions_same_spin(p, q, r)-= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* current.epsilon_tilde[p]
	* occupation_numbers[p];
dHdl.interactions_same_spin(p, q, r)-= 2.000000 * sign(current.epsilon_tilde[p+r] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(p+r, q-r, -r)
	* current.epsilon_tilde[q]
	* occupation_numbers[q];
for (momentum_iterator<L> s = momentum_iterator<L>::begin(); s != momentum_iterator<L>::end(); ++s) {
double ns_value{};
double one_value{};
one_value += 2.000000 * current.interactions_same_spin(p, q, s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p-s] - current.epsilon_tilde[-q+s]) 
	* current.interactions_same_spin(-p-r, -q+r, r-s);
one_value -= 2.000000 * current.interactions_same_spin(-s, p+q+s, p+r+s)
	* sign(current.epsilon_tilde[-s] + current.epsilon_tilde[p+q+s] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(-s, p+q+s, p+s);
one_value += 2.000000 * current.interactions_same_spin(-s, p+q+s, p+r+s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[s] - current.epsilon_tilde[-p-q-s]) 
	* current.interactions_same_spin(-p, -q, p+s);
one_value -= 2.000000 * current.interactions_same_spin(p, q, s)
	* sign(current.epsilon_tilde[p+s] + current.epsilon_tilde[q-s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p+s, q-s, r-s);
ns_value -= 2.000000 * current.interactions_same_spin(p, q, q-s)
	* sign(current.epsilon_tilde[p+r] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p+q-s] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(p+r, q-r, q-r-s);
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r);
ns_value += 2.000000 * current.interactions_same_spin(-p, -q, p-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+q-s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(s, p+q-s, p+r-s);
ns_value -= 2.000000 * current.interactions_same_spin(-p, -q, p-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p-q+s]) 
	* current.interactions_same_spin(-p-r, -q+r, p+r-s);
ns_value += 2.000000 * current.interactions_same_spin(p, q, q-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+q-s] - current.epsilon_tilde[q-r] - current.epsilon_tilde[p+r]) 
	* current.interactions_same_spin(s, p+q-s, q-r-s);
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_same_spin(-p, -r-s, p-s);
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, p+r-s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_same_spin(-p, p-q+r-s, p-s);
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, q-r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(p, r+s, r);
ns_value -= 2.000000 * current.interactions_differing_spin(s, q-r, GammaPoint<L>)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p-r, -q+r, r);
ns_value += 4.000000 * current.interactions_differing_spin(-q, -s, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[r+s] - current.epsilon_tilde[p]) 
	* current.interactions_differing_spin(s, p+r, r);
ns_value -= 2.000000 * current.interactions_differing_spin(s, p+r, GammaPoint<L>)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p-r, -q+r, r);
ns_value -= 4.000000 * current.interactions_differing_spin(q, s, GammaPoint<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r);
ns_value -= 4.000000 * current.interactions_differing_spin(p, r-s, r)
	* sign(current.epsilon_tilde[q-r] + current.epsilon_tilde[r-s] - current.epsilon_tilde[q] - current.epsilon_tilde[-s]) 
	* current.interactions_differing_spin(q-r, r-s, r);
ns_value += 4.000000 * current.interactions_differing_spin(p, r-s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-r+s]) 
	* current.interactions_differing_spin(-q, s, r);
ns_value -= 4.000000 * current.interactions_differing_spin(-q, -s, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, r+s, r);
ns_value += 4.000000 * current.interactions_differing_spin(q, s, GammaPoint<L>)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p-r, -q+r, r);
ns_value += 2.000000 * current.interactions_differing_spin(s, p+r, GammaPoint<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r);
ns_value += 2.000000 * current.interactions_differing_spin(s, q-r, GammaPoint<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r);
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, p+r-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[-p+q-r+s]) 
	* current.interactions_same_spin(s, q-r, p-s);
ns_value += 4.000000 * current.interactions_same_spin(-p, p-q+r-s, p-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p-r, p-q+r-s, p+r-s);
ns_value += 4.000000 * current.interactions_same_spin(-p, -r+s, p+s)
	* sign(current.epsilon_tilde[q-r] + current.epsilon_tilde[r-s] - current.epsilon_tilde[q] - current.epsilon_tilde[-s]) 
	* current.interactions_same_spin(q-r, r-s, r);
ns_value += 4.000000 * current.interactions_same_spin(p, r-s, r)
	* sign(current.epsilon_tilde[-q+r] + current.epsilon_tilde[-r+s] - current.epsilon_tilde[s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q+r, -r+s, q-r+s);
ns_value -= 2.000000 * current.interactions_same_spin(-s, p+q+s, p+r+s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[s] - current.epsilon_tilde[-p-q-s]) 
	* current.interactions_same_spin(-p, -q, p+s);
ns_value -= 2.000000 * current.interactions_same_spin(s, q-r, q-r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r);
ns_value -= 2.000000 * current.interactions_same_spin(s, p+r, p+r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r);
ns_value -= 2.000000 * current.interactions_same_spin(-s, p+q+s, q-r+s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+q+s] - current.epsilon_tilde[-s]) 
	* current.interactions_same_spin(p, q, q+s);
ns_value -= 4.000000 * current.interactions_same_spin(-p, -r+s, p+s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-r+s]) 
	* current.interactions_same_spin(-q, s, r);
ns_value -= 4.000000 * current.interactions_same_spin(-p, p-q+r-s, p-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-s] - current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_same_spin(-q, -s, p+r-s);
ns_value -= 4.000000 * current.interactions_same_spin(p, r-s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[-r+s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_same_spin(-q, s, q-r+s);
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, q-r-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[r+s] - current.epsilon_tilde[p]) 
	* current.interactions_same_spin(s, p+r, r);
ns_value += 2.000000 * current.interactions_same_spin(-s, p+q+s, q-r+s)
	* sign(current.epsilon_tilde[-s] + current.epsilon_tilde[p+q+s] - current.epsilon_tilde[q] - current.epsilon_tilde[p]) 
	* current.interactions_same_spin(-s, p+q+s, q+s);
ns_value += 2.000000 * current.interactions_same_spin(s, p+r, p+r-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p-r, -q+r, r);
ns_value += 2.000000 * current.interactions_same_spin(s, q-r, q-r-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p-r, -q+r, r);
ns_value += 2.000000 * current.interactions_same_spin(-s, p+q+s, p+r+s)
	* sign(current.epsilon_tilde[-s] + current.epsilon_tilde[p+q+s] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(-s, p+q+s, p+s);
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[p] - current.epsilon_tilde[r+s]) 
	* current.interactions_same_spin(s, p+r, p-s);
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p-r, -q+r, r);
ns_value -= 4.000000 * current.interactions_differing_spin(q, s, GammaPoint<L>)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p-r, -q+r, r)
	* occupation_numbers[q];
ns_value -= 4.000000 * current.interactions_differing_spin(p, s, GammaPoint<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* occupation_numbers[p];
ns_value -= 4.000000 * current.interactions_differing_spin(p, s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r-s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-s]) 
	* current.interactions_differing_spin(-q, r-s, r)
	* occupation_numbers[r-s];
ns_value += 4.000000 * current.interactions_differing_spin(p, s, GammaPoint<L>)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p-r, -q+r, r)
	* occupation_numbers[p];
ns_value += 4.000000 * current.interactions_differing_spin(p, s, r)
	* sign(current.epsilon_tilde[-s] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[r-s] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-s, -q+r, r)
	* occupation_numbers[r-s];
ns_value -= 4.000000 * current.interactions_differing_spin(-q, -s, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[r+s] - current.epsilon_tilde[p]) 
	* current.interactions_differing_spin(s, p+r, r)
	* occupation_numbers[r+s];
ns_value += 4.000000 * current.interactions_differing_spin(q, s, GammaPoint<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* occupation_numbers[q];
ns_value += 4.000000 * current.interactions_differing_spin(-q, -s, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, r+s, r)
	* occupation_numbers[r+s];
ns_value -= 4.000000 * current.interactions_same_spin(-p, -s, p-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p-r, -q+r, r)
	* occupation_numbers[p];
ns_value -= 4.000000 * current.interactions_same_spin(-p, p-q+r-s, p-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p-r, p-q+r-s, p+r-s)
	* occupation_numbers[-p+q-r+s];
ns_value += 4.000000 * current.interactions_same_spin(-p, -s, p+r-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r-s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-s]) 
	* current.interactions_same_spin(-q, r-s, r)
	* occupation_numbers[r-s];
ns_value -= 4.000000 * current.interactions_same_spin(p, s, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[q] - current.epsilon_tilde[-r+s]) 
	* current.interactions_same_spin(s, q-r, q-s)
	* occupation_numbers[r-s];
ns_value -= 4.000000 * current.interactions_same_spin(-p, -s, p+r-s)
	* sign(current.epsilon_tilde[-s] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[r-s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-s, -q+r, r)
	* occupation_numbers[r-s];
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, p+r-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[-p+q-r+s]) 
	* current.interactions_same_spin(s, q-r, p-s)
	* occupation_numbers[p-q+r-s];
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, q-r-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[r+s] - current.epsilon_tilde[p]) 
	* current.interactions_same_spin(s, p+r, r)
	* occupation_numbers[r+s];
ns_value += 4.000000 * current.interactions_same_spin(-p, -s, p-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* occupation_numbers[p];
ns_value += 4.000000 * current.interactions_same_spin(-p, p-q+r-s, p-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-s] - current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_same_spin(-q, -s, p+r-s)
	* occupation_numbers[-p+q-r+s];
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p-r, -q+r, r)
	* occupation_numbers[q];
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, p+r-s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_same_spin(-p, p-q+r-s, p-s)
	* occupation_numbers[p-q+r-s];
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, q-r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(p, r+s, r)
	* occupation_numbers[r+s];
ns_value += 4.000000 * current.interactions_same_spin(p, s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_same_spin(-q, r-s, q-s)
	* occupation_numbers[r-s];
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_same_spin(-p, -r-s, p-s)
	* occupation_numbers[r+s];
ns_value -= 4.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* occupation_numbers[q];
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[p] - current.epsilon_tilde[r+s]) 
	* current.interactions_same_spin(s, p+r, p-s)
	* occupation_numbers[r+s];
dHdl.interactions_same_spin(p, q, r) += one_value + occupation_numbers[s] * ns_value;
} // s-loop
} // r-loop
} // q-loop
} // p-loop
dHdl.interactions_same_spin.symmetrize();
dHdl.interactions_differing_spin.symmetrize();
}
} // namespace NickelCUT::flow