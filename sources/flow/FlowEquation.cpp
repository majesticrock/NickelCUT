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
one_value += 8.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-q] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p-r, -q+r, p-q+r);
one_value += 8.000000 * current.interactions_differing_spin(-r, -p+q+r, q+r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[q] - current.epsilon_tilde[-p+q+r] - current.epsilon_tilde[-r]) 
	* current.interactions_differing_spin(-p, q, q+r);
one_value -= 8.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[-p-r] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-p] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p-r, -q+r, r);
one_value -= 8.000000 * current.interactions_same_spin(-r, -p+q+r, q+r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[q] - current.epsilon_tilde[-r] - current.epsilon_tilde[-p+q+r]) 
	* current.interactions_same_spin(-p, q, p-r);
nr_value += 2.000000 * current.interactions_same_spin(-q, -r, q-r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[-r] - current.epsilon_tilde[p] - current.epsilon_tilde[-r]) 
	* current.interactions_differing_spin(p, -r, Gamma<L>);
nr_value -= 8.000000 * current.interactions_differing_spin(-p, p-q-r, p-q)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[-p+q+r]) 
	* current.interactions_differing_spin(q, r, p-q);
nr_value += 8.000000 * current.interactions_differing_spin(-p, -q, p-r)
	* sign(current.epsilon_tilde[r] + current.epsilon_tilde[p+q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_differing_spin(r, p+q-r, p-r);
nr_value += 2.000000 * current.interactions_differing_spin(p, q, Gamma<L>)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r] - current.epsilon_tilde[-r] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -r, p-r);
nr_value += 4.000000 * current.interactions_differing_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[-q] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -q, p-q);
nr_value -= 8.000000 * current.interactions_differing_spin(-q, -p+q+r, q+r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[r] - current.epsilon_tilde[-p+q+r] - current.epsilon_tilde[-q]) 
	* current.interactions_differing_spin(-p, r, q+r);
nr_value -= 8.000000 * current.interactions_differing_spin(q, p-q-r, p-q)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[r] - current.epsilon_tilde[-q] - current.epsilon_tilde[-p+q+r]) 
	* current.interactions_differing_spin(-p, r, p-q);
nr_value -= 2.000000 * current.interactions_differing_spin(p, q, Gamma<L>)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[r] - current.epsilon_tilde[r] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, r, p+r);
nr_value -= 2.000000 * current.interactions_differing_spin(p, q, Gamma<L>)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r] - current.epsilon_tilde[r] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, r, q+r);
nr_value -= 6.000000 * current.interactions_same_spin(-q, -r, q-r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(p, r, Gamma<L>);
nr_value -= 4.000000 * current.interactions_differing_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r] - current.epsilon_tilde[-r] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -r, p-r);
nr_value -= 4.000000 * current.interactions_differing_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_differing_spin(p, q, Gamma<L>);
nr_value += 4.000000 * current.interactions_differing_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(p, r, Gamma<L>);
nr_value += 8.000000 * current.interactions_differing_spin(q, r, p-q)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q-r] - current.epsilon_tilde[-q] - current.epsilon_tilde[-r]) 
	* current.interactions_differing_spin(-p, p-q-r, p-q);
nr_value += 4.000000 * current.interactions_same_spin(-q, -r, q-r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r] - current.epsilon_tilde[-r] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -r, p-r);
nr_value += 4.000000 * current.interactions_same_spin(-q, -r, q-r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_differing_spin(p, q, Gamma<L>);
nr_value += 8.000000 * current.interactions_same_spin(-q, -p+q+r, q+r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[r] - current.epsilon_tilde[-q] - current.epsilon_tilde[-p+q+r]) 
	* current.interactions_same_spin(-p, r, p-q);
nr_value += 8.000000 * current.interactions_same_spin(q, p-q-r, p-q)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[r] - current.epsilon_tilde[-p+q+r] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p, r, q+r);
nr_value += 8.000000 * current.interactions_same_spin(-p, p-q-r, p-q)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-r] - current.epsilon_tilde[p-q-r] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-q, -r, p-r);
nr_value -= 8.000000 * current.interactions_same_spin(q, r, p-q)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q-r] - current.epsilon_tilde[-r] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-p, p-q-r, p-r);
nr_value -= 8.000000 * current.interactions_same_spin(p, q, q-r)
	* sign(current.epsilon_tilde[r] + current.epsilon_tilde[p+q-r] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(r, p+q-r, p-r);
nr_value -= 8.000000 * current.interactions_same_spin(-p, -q, p-r)
	* sign(current.epsilon_tilde[r] + current.epsilon_tilde[p+q-r] - current.epsilon_tilde[q] - current.epsilon_tilde[p]) 
	* current.interactions_same_spin(r, p+q-r, q-r);
nr_value -= 4.000000 * current.interactions_same_spin(-q, -r, q-r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[-q] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -q, p-q);
nr_value += 8.000000 * current.interactions_differing_spin(p, q, q-r)
	* sign(current.epsilon_tilde[r] + current.epsilon_tilde[p+q-r] - current.epsilon_tilde[q] - current.epsilon_tilde[p]) 
	* current.interactions_differing_spin(r, p+q-r, q-r);
nr_value += 2.000000 * current.interactions_differing_spin(p, q, Gamma<L>)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-r] - current.epsilon_tilde[-r] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, -r, q-r);
nr_value += 4.000000 * current.interactions_differing_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_differing_spin(p, q, Gamma<L>);
nr_value -= 4.000000 * current.interactions_differing_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(p, r, Gamma<L>);
nr_value -= 2.000000 * current.interactions_differing_spin(p, q, Gamma<L>)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-r] - current.epsilon_tilde[-r] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, -r, q-r);
nr_value += 4.000000 * current.interactions_differing_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r] - current.epsilon_tilde[-r] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -r, p-r);
nr_value += 4.000000 * current.interactions_same_spin(-q, -r, q-r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(p, r, Gamma<L>);
nr_value -= 4.000000 * current.interactions_differing_spin(p, q, Gamma<L>)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r] - current.epsilon_tilde[-r] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -r, p-r);
nr_value -= 4.000000 * current.interactions_differing_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[-q] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -q, p-q);
nr_value += 2.000000 * current.interactions_differing_spin(p, q, Gamma<L>)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r] - current.epsilon_tilde[r] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, r, q+r);
nr_value -= 4.000000 * current.interactions_same_spin(-q, -r, q-r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p] - current.epsilon_tilde[q]) 
	* current.interactions_differing_spin(p, q, Gamma<L>);
nr_value += 4.000000 * current.interactions_differing_spin(p, q, Gamma<L>)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[r] - current.epsilon_tilde[r] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, r, p+r);
nr_value -= 4.000000 * current.interactions_same_spin(-q, -r, q-r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r] - current.epsilon_tilde[-r] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -r, p-r);
nr_value += 4.000000 * current.interactions_same_spin(-q, -r, q-r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[-q] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -q, p-q);
nr_value -= 2.000000 * current.interactions_differing_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[s] - current.epsilon_tilde[s] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, s, p+s);
nr_value += 2.000000 * current.interactions_differing_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -s, p-s);
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
dHdl.interactions_differing_spin(p, q, r)-= 4.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* current.epsilon_tilde[q];
dHdl.interactions_differing_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* current.epsilon_tilde[q-r];
dHdl.interactions_differing_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* current.epsilon_tilde[p+r];
dHdl.interactions_differing_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[p+r] - current.epsilon_tilde[p+r] - current.epsilon_tilde[p]) 
	* current.interactions_same_spin(p, p+r, r)
	* current.epsilon_tilde[q]
	* occupation_numbers[p+r];
dHdl.interactions_differing_spin(p, q, r)-= 2.000000 * sign(current.epsilon_tilde[q] + current.epsilon_tilde[q+r] - current.epsilon_tilde[q+r] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(q, q+r, r)
	* current.epsilon_tilde[p]
	* occupation_numbers[q+r];
dHdl.interactions_differing_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[p-r] - current.epsilon_tilde[p-r] - current.epsilon_tilde[p]) 
	* current.interactions_same_spin(p, p-r, -r)
	* current.epsilon_tilde[q]
	* occupation_numbers[p-r];
dHdl.interactions_differing_spin(p, q, r)-= 4.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* current.epsilon_tilde[p]
	* occupation_numbers[p];
dHdl.interactions_differing_spin(p, q, r)-= 2.000000 * sign(current.epsilon_tilde[q] + current.epsilon_tilde[q-r] - current.epsilon_tilde[q-r] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(q, q-r, -r)
	* current.epsilon_tilde[p]
	* occupation_numbers[q-r];
dHdl.interactions_differing_spin(p, q, r)+= 4.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* current.epsilon_tilde[q]
	* occupation_numbers[q];
dHdl.interactions_differing_spin(p, q, r)-= 4.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(p, r, Gamma<L>)
	* current.epsilon_tilde[q]
	* occupation_numbers[r];
dHdl.interactions_differing_spin(p, q, r)+= 4.000000 * sign(current.epsilon_tilde[q] + current.epsilon_tilde[r] - current.epsilon_tilde[q] - current.epsilon_tilde[r]) 
	* current.interactions_same_spin(q, r, Gamma<L>)
	* current.epsilon_tilde[p]
	* occupation_numbers[r];
for (momentum_iterator<L> s = momentum_iterator<L>::begin(); s != momentum_iterator<L>::end(); ++s) {
double ns_value{};
double one_value{};
one_value += 4.000000 * current.interactions_differing_spin(-s, p+q+s, p+r+s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[s] - current.epsilon_tilde[-p-q-s]) 
	* current.interactions_differing_spin(-p, -q, p+s);
one_value -= 4.000000 * current.interactions_differing_spin(p, q, s)
	* sign(current.epsilon_tilde[p+s] + current.epsilon_tilde[q-s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p+s, q-s, r-s);
ns_value -= 8.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, -s, q-s);
ns_value += 4.000000 * current.interactions_same_spin(s, q-r, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r);
ns_value -= 4.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[s] - current.epsilon_tilde[p+r]) 
	* current.interactions_differing_spin(s, p+r, Gamma<L>);
ns_value -= 8.000000 * current.interactions_same_spin(-q, -s, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, r+s, r);
ns_value += 8.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[s] - current.epsilon_tilde[q] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(q, s, Gamma<L>);
ns_value -= 8.000000 * current.interactions_same_spin(q, s, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r);
ns_value += 8.000000 * current.interactions_differing_spin(p, r-s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-r+s]) 
	* current.interactions_same_spin(-q, s, r);
ns_value += 4.000000 * current.interactions_differing_spin(s, p+r, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r);
ns_value -= 4.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[s] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(s, q-r, Gamma<L>);
ns_value -= 8.000000 * current.interactions_same_spin(-p, -r+s, p+s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-r+s]) 
	* current.interactions_differing_spin(-q, s, r);
ns_value -= 4.000000 * current.interactions_differing_spin(-s, p+q+s, q-r+s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+q+s] - current.epsilon_tilde[-s]) 
	* current.interactions_differing_spin(p, q, q+s);
ns_value -= 8.000000 * current.interactions_differing_spin(-p, p-q+r-s, p-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-s] - current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_differing_spin(-q, -s, p+r-s);
ns_value -= 4.000000 * current.interactions_same_spin(s, q-r, q-r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r);
ns_value += 8.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r);
ns_value += 8.000000 * current.interactions_differing_spin(-q, -s, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_same_spin(-p, -r-s, p-s);
ns_value += 4.000000 * current.interactions_differing_spin(p, q, q-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+q-s] - current.epsilon_tilde[q-r] - current.epsilon_tilde[p+r]) 
	* current.interactions_differing_spin(s, p+q-s, q-r-s);
ns_value += 4.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(s, p+r, p+r-s);
ns_value += 4.000000 * current.interactions_differing_spin(-p, -q, p-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+q-s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(s, p+q-s, p+r-s);
ns_value -= 8.000000 * current.interactions_differing_spin(p, r-s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[-r+s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_same_spin(-q, s, q-r+s);
ns_value += 4.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[q-r] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(s, q-r, q-r-s);
ns_value -= 4.000000 * current.interactions_differing_spin(-s, p+q+s, p+r+s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[s] - current.epsilon_tilde[-p-q-s]) 
	* current.interactions_differing_spin(-p, -q, p+s);
ns_value += 8.000000 * current.interactions_differing_spin(-q, -s, p+r-s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_differing_spin(-p, p-q+r-s, p-s);
ns_value -= 4.000000 * current.interactions_same_spin(s, p+r, p+r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r);
ns_value += 8.000000 * current.interactions_same_spin(-q, -s, q-r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, r+s, r);
ns_value += 8.000000 * current.interactions_same_spin(-q, -s, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, r+s, r)
	* occupation_numbers[r+s];
ns_value += 8.000000 * current.interactions_same_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -s, p-s)
	* occupation_numbers[r];
ns_value -= 2.000000 * current.interactions_differing_spin(p, s, Gamma<L>)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[q+r] - current.epsilon_tilde[q+r] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(q, q+r, r)
	* occupation_numbers[q+r];
ns_value -= 2.000000 * current.interactions_differing_spin(p, r, Gamma<L>)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, -s, q-s)
	* occupation_numbers[r];
ns_value -= 8.000000 * current.interactions_same_spin(-p, -s, p-s)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[r] - current.epsilon_tilde[q] - current.epsilon_tilde[r]) 
	* current.interactions_same_spin(q, r, Gamma<L>)
	* occupation_numbers[r];
ns_value -= 2.000000 * current.interactions_differing_spin(p, r, Gamma<L>)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, s, q+s)
	* occupation_numbers[r];
ns_value += 2.000000 * current.interactions_same_spin(-q, -q+r, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[s] - current.epsilon_tilde[p] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, s, Gamma<L>)
	* occupation_numbers[q-r];
ns_value -= 8.000000 * current.interactions_differing_spin(p, s, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* occupation_numbers[p];
ns_value += 8.000000 * current.interactions_same_spin(q, s, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* occupation_numbers[q];
ns_value += 8.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[s] - current.epsilon_tilde[p] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, s, Gamma<L>)
	* occupation_numbers[p];
ns_value += 2.000000 * current.interactions_same_spin(q, q+r, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[s] - current.epsilon_tilde[p] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, s, Gamma<L>)
	* occupation_numbers[q+r];
ns_value += 2.000000 * current.interactions_differing_spin(p, s, Gamma<L>)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[r] - current.epsilon_tilde[q] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(q, r, Gamma<L>)
	* occupation_numbers[r];
ns_value += 2.000000 * current.interactions_differing_spin(p, s, Gamma<L>)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[r] - current.epsilon_tilde[q] - current.epsilon_tilde[r]) 
	* current.interactions_same_spin(q, r, Gamma<L>)
	* occupation_numbers[r];
ns_value -= 2.000000 * current.interactions_same_spin(q, s, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(p, r, Gamma<L>)
	* occupation_numbers[r];
ns_value -= 2.000000 * current.interactions_differing_spin(q, s, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(p, r, Gamma<L>)
	* occupation_numbers[r];
ns_value += 2.000000 * current.interactions_differing_spin(p, r, Gamma<L>)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[s] - current.epsilon_tilde[q] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(q, s, Gamma<L>)
	* occupation_numbers[r];
ns_value -= 2.000000 * current.interactions_same_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[s] - current.epsilon_tilde[p] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, s, Gamma<L>)
	* occupation_numbers[r];
ns_value += 2.000000 * current.interactions_differing_spin(p, r, Gamma<L>)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[s] - current.epsilon_tilde[q] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(q, s, Gamma<L>)
	* occupation_numbers[r];
ns_value -= 2.000000 * current.interactions_differing_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[s] - current.epsilon_tilde[p] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, s, Gamma<L>)
	* occupation_numbers[r];
ns_value -= 8.000000 * current.interactions_same_spin(-q, -s, q-r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, r+s, r)
	* occupation_numbers[r+s];
ns_value += 8.000000 * current.interactions_differing_spin(-p, -s, q-r-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-p+q-r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_differing_spin(-q, -p+q-r-s, q-s)
	* occupation_numbers[p-q+r+s];
ns_value += 4.000000 * current.interactions_same_spin(-p, -s, p-s)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[q+r] - current.epsilon_tilde[q+r] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(q, q+r, r)
	* occupation_numbers[q+r];
ns_value += 8.000000 * current.interactions_differing_spin(p, s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_same_spin(-q, r-s, q-s)
	* occupation_numbers[r-s];
ns_value += 8.000000 * current.interactions_same_spin(-p, -s, p-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* occupation_numbers[p];
ns_value -= 2.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[q-r] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(s, q-r, q-r-s)
	* occupation_numbers[q-r];
ns_value += 4.000000 * current.interactions_same_spin(-p, -s, p-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, -q+r, r)
	* occupation_numbers[q-r];
ns_value += 8.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, -s, q-s)
	* occupation_numbers[q];
ns_value += 8.000000 * current.interactions_same_spin(-p, -s, p+r-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r-s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-s]) 
	* current.interactions_differing_spin(-q, r-s, r)
	* occupation_numbers[r-s];
ns_value += 2.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[-s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[q-r] - current.epsilon_tilde[-s]) 
	* current.interactions_same_spin(-s, q-r, q-r+s)
	* occupation_numbers[q-r];
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(p, r, Gamma<L>)
	* occupation_numbers[r];
ns_value -= 8.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(p, q, r)
	* occupation_numbers[q];
ns_value -= 8.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -s, p-s)
	* occupation_numbers[p];
ns_value -= 4.000000 * current.interactions_same_spin(q, q+r, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -s, p-s)
	* occupation_numbers[q+r];
ns_value -= 4.000000 * current.interactions_same_spin(-q, -q+r, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -s, p-s)
	* occupation_numbers[q-r];
ns_value -= 8.000000 * current.interactions_differing_spin(-q, -s, p+r-s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_differing_spin(-p, p-q+r-s, p-s)
	* occupation_numbers[p-q+r-s];
ns_value -= 8.000000 * current.interactions_differing_spin(-q, -s, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_same_spin(-p, -r-s, p-s)
	* occupation_numbers[r+s];
ns_value -= 8.000000 * current.interactions_differing_spin(p, q, r)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[s] - current.epsilon_tilde[q] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(q, s, Gamma<L>)
	* occupation_numbers[q];
ns_value -= 8.000000 * current.interactions_differing_spin(p, s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r-s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-s]) 
	* current.interactions_same_spin(-q, r-s, r)
	* occupation_numbers[r-s];
ns_value -= 2.000000 * current.interactions_differing_spin(p, s, Gamma<L>)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, -q+r, r)
	* occupation_numbers[q-r];
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
dHdl.interactions_same_spin(p, q, r)-= 4.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* current.epsilon_tilde[q];
dHdl.interactions_same_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* current.epsilon_tilde[q-r];
dHdl.interactions_same_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* current.epsilon_tilde[p+r];
dHdl.interactions_same_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[p+r] - current.epsilon_tilde[p+r] - current.epsilon_tilde[p]) 
	* current.interactions_same_spin(p, p+r, r)
	* current.epsilon_tilde[q]
	* occupation_numbers[p+r];
dHdl.interactions_same_spin(p, q, r)-= 2.000000 * sign(current.epsilon_tilde[q] + current.epsilon_tilde[q+r] - current.epsilon_tilde[q+r] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(q, q+r, r)
	* current.epsilon_tilde[p]
	* occupation_numbers[q+r];
dHdl.interactions_same_spin(p, q, r)+= 2.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[p-r] - current.epsilon_tilde[p-r] - current.epsilon_tilde[p]) 
	* current.interactions_same_spin(p, p-r, -r)
	* current.epsilon_tilde[q]
	* occupation_numbers[p-r];
dHdl.interactions_same_spin(p, q, r)-= 4.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* current.epsilon_tilde[p]
	* occupation_numbers[p];
dHdl.interactions_same_spin(p, q, r)-= 2.000000 * sign(current.epsilon_tilde[q] + current.epsilon_tilde[q-r] - current.epsilon_tilde[q-r] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(q, q-r, -r)
	* current.epsilon_tilde[p]
	* occupation_numbers[q-r];
dHdl.interactions_same_spin(p, q, r)+= 4.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* current.epsilon_tilde[q]
	* occupation_numbers[q];
dHdl.interactions_same_spin(p, q, r)-= 4.000000 * sign(current.epsilon_tilde[p] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(p, r, Gamma<L>)
	* current.epsilon_tilde[q]
	* occupation_numbers[r];
dHdl.interactions_same_spin(p, q, r)+= 4.000000 * sign(current.epsilon_tilde[q] + current.epsilon_tilde[r] - current.epsilon_tilde[q] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(q, r, Gamma<L>)
	* current.epsilon_tilde[p]
	* occupation_numbers[r];
for (momentum_iterator<L> s = momentum_iterator<L>::begin(); s != momentum_iterator<L>::end(); ++s) {
double ns_value{};
double one_value{};
one_value += 4.000000 * current.interactions_same_spin(-s, p+q+s, p+r+s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[s] - current.epsilon_tilde[-p-q-s]) 
	* current.interactions_same_spin(-p, -q, p+s);
one_value -= 4.000000 * current.interactions_same_spin(p, q, s)
	* sign(current.epsilon_tilde[p+s] + current.epsilon_tilde[q-s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p+s, q-s, r-s);
ns_value -= 8.000000 * current.interactions_same_spin(-p, -r+s, p+s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-r+s]) 
	* current.interactions_same_spin(-q, s, r);
ns_value += 4.000000 * current.interactions_differing_spin(s, q-r, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r);
ns_value -= 4.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[s] - current.epsilon_tilde[p+r]) 
	* current.interactions_differing_spin(s, p+r, Gamma<L>);
ns_value -= 8.000000 * current.interactions_differing_spin(-q, -s, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, r+s, r);
ns_value += 8.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[s] - current.epsilon_tilde[q] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(q, s, Gamma<L>);
ns_value -= 8.000000 * current.interactions_differing_spin(q, s, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r);
ns_value += 8.000000 * current.interactions_differing_spin(p, r-s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-r+s]) 
	* current.interactions_differing_spin(-q, s, r);
ns_value += 4.000000 * current.interactions_differing_spin(s, p+r, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r);
ns_value -= 4.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[s] - current.epsilon_tilde[q-r]) 
	* current.interactions_differing_spin(s, q-r, Gamma<L>);
ns_value -= 8.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, -s, q-s);
ns_value -= 4.000000 * current.interactions_same_spin(-s, p+q+s, q-r+s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+q+s] - current.epsilon_tilde[-s]) 
	* current.interactions_same_spin(p, q, q+s);
ns_value -= 8.000000 * current.interactions_same_spin(-p, p-q+r-s, p-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-s] - current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_same_spin(-q, -s, p+r-s);
ns_value -= 4.000000 * current.interactions_same_spin(s, q-r, q-r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r);
ns_value -= 4.000000 * current.interactions_same_spin(s, p+r, p+r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r);
ns_value += 8.000000 * current.interactions_same_spin(-q, -s, q-r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(p, r+s, r);
ns_value += 8.000000 * current.interactions_same_spin(-q, -s, p+r-s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_same_spin(-p, p-q+r-s, p-s);
ns_value -= 4.000000 * current.interactions_same_spin(-s, p+q+s, p+r+s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-q] - current.epsilon_tilde[s] - current.epsilon_tilde[-p-q-s]) 
	* current.interactions_same_spin(-p, -q, p+s);
ns_value += 4.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[q-r] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(s, q-r, q-r-s);
ns_value -= 8.000000 * current.interactions_same_spin(p, r-s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[-r+s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_same_spin(-q, s, q-r+s);
ns_value += 4.000000 * current.interactions_same_spin(-p, -q, p-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+q-s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(s, p+q-s, p+r-s);
ns_value += 4.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+r] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(s, p+r, p+r-s);
ns_value += 4.000000 * current.interactions_same_spin(p, q, q-s)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[p+q-s] - current.epsilon_tilde[q-r] - current.epsilon_tilde[p+r]) 
	* current.interactions_same_spin(s, p+q-s, q-r-s);
ns_value += 8.000000 * current.interactions_same_spin(-q, -s, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_same_spin(-p, -r-s, p-s);
ns_value += 8.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r);
ns_value -= 2.000000 * current.interactions_differing_spin(p, s, Gamma<L>)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, -q+r, r)
	* occupation_numbers[q-r];
ns_value += 4.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(p, r, Gamma<L>)
	* occupation_numbers[r];
ns_value -= 2.000000 * current.interactions_differing_spin(p, s, Gamma<L>)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[q+r] - current.epsilon_tilde[q+r] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(q, q+r, r)
	* occupation_numbers[q+r];
ns_value += 8.000000 * current.interactions_differing_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -s, p-s)
	* occupation_numbers[r];
ns_value += 8.000000 * current.interactions_differing_spin(-q, -s, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, r+s, r)
	* occupation_numbers[r+s];
ns_value -= 2.000000 * current.interactions_differing_spin(p, r, Gamma<L>)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, -s, q-s)
	* occupation_numbers[r];
ns_value -= 8.000000 * current.interactions_same_spin(-p, -s, p-s)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[r] - current.epsilon_tilde[q] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(q, r, Gamma<L>)
	* occupation_numbers[r];
ns_value -= 2.000000 * current.interactions_differing_spin(p, r, Gamma<L>)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[s] - current.epsilon_tilde[s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, s, q+s)
	* occupation_numbers[r];
ns_value += 2.000000 * current.interactions_same_spin(-q, -q+r, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[s] - current.epsilon_tilde[p] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, s, Gamma<L>)
	* occupation_numbers[q-r];
ns_value -= 8.000000 * current.interactions_differing_spin(p, s, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* occupation_numbers[p];
ns_value += 8.000000 * current.interactions_differing_spin(q, s, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* occupation_numbers[q];
ns_value += 8.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[s] - current.epsilon_tilde[p] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, s, Gamma<L>)
	* occupation_numbers[p];
ns_value += 2.000000 * current.interactions_same_spin(q, q+r, r)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[s] - current.epsilon_tilde[p] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, s, Gamma<L>)
	* occupation_numbers[q+r];
ns_value += 4.000000 * current.interactions_differing_spin(p, s, Gamma<L>)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[r] - current.epsilon_tilde[q] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(q, r, Gamma<L>)
	* occupation_numbers[r];
ns_value -= 4.000000 * current.interactions_differing_spin(q, s, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r] - current.epsilon_tilde[p] - current.epsilon_tilde[r]) 
	* current.interactions_differing_spin(p, r, Gamma<L>)
	* occupation_numbers[r];
ns_value += 4.000000 * current.interactions_differing_spin(p, r, Gamma<L>)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[s] - current.epsilon_tilde[q] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(q, s, Gamma<L>)
	* occupation_numbers[r];
ns_value -= 4.000000 * current.interactions_differing_spin(q, r, Gamma<L>)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[s] - current.epsilon_tilde[p] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(p, s, Gamma<L>)
	* occupation_numbers[r];
ns_value += 2.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[-s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[q-r] - current.epsilon_tilde[-s]) 
	* current.interactions_same_spin(-s, q-r, q-r+s)
	* occupation_numbers[q-r];
ns_value += 8.000000 * current.interactions_same_spin(-p, -s, q-r-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-p+q-r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_same_spin(-q, -p+q-r-s, q-s)
	* occupation_numbers[p-q+r+s];
ns_value += 4.000000 * current.interactions_same_spin(-p, -s, p-s)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[q+r] - current.epsilon_tilde[q+r] - current.epsilon_tilde[q]) 
	* current.interactions_same_spin(q, q+r, r)
	* occupation_numbers[q+r];
ns_value += 8.000000 * current.interactions_same_spin(p, s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_same_spin(-q, r-s, q-s)
	* occupation_numbers[r-s];
ns_value += 8.000000 * current.interactions_same_spin(-p, -s, p-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* occupation_numbers[p];
ns_value -= 2.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[s] + current.epsilon_tilde[q-r] - current.epsilon_tilde[q-r] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(s, q-r, q-r-s)
	* occupation_numbers[q-r];
ns_value += 4.000000 * current.interactions_same_spin(-p, -s, p-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-q+r] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, -q+r, r)
	* occupation_numbers[q-r];
ns_value += 8.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q]) 
	* current.interactions_same_spin(-q, -s, q-s)
	* occupation_numbers[q];
ns_value += 8.000000 * current.interactions_same_spin(-p, -s, p+r-s)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r-s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-s]) 
	* current.interactions_same_spin(-q, r-s, r)
	* occupation_numbers[r-s];
ns_value -= 8.000000 * current.interactions_differing_spin(p, s, r)
	* sign(current.epsilon_tilde[-q] + current.epsilon_tilde[r-s] - current.epsilon_tilde[-q+r] - current.epsilon_tilde[-s]) 
	* current.interactions_differing_spin(-q, r-s, r)
	* occupation_numbers[r-s];
ns_value -= 8.000000 * current.interactions_same_spin(-q, -s, q-r-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[r+s] - current.epsilon_tilde[p+r] - current.epsilon_tilde[s]) 
	* current.interactions_same_spin(p, r+s, r)
	* occupation_numbers[r+s];
ns_value -= 8.000000 * current.interactions_same_spin(-q, -s, q-s)
	* sign(current.epsilon_tilde[p] + current.epsilon_tilde[q] - current.epsilon_tilde[p+r] - current.epsilon_tilde[q-r]) 
	* current.interactions_same_spin(p, q, r)
	* occupation_numbers[q];
ns_value -= 8.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -s, p-s)
	* occupation_numbers[p];
ns_value -= 4.000000 * current.interactions_same_spin(q, q+r, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -s, p-s)
	* occupation_numbers[q+r];
ns_value -= 4.000000 * current.interactions_same_spin(-q, -q+r, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p]) 
	* current.interactions_same_spin(-p, -s, p-s)
	* occupation_numbers[q-r];
ns_value -= 8.000000 * current.interactions_same_spin(-q, -s, p+r-s)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[p-q+r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-q+r]) 
	* current.interactions_same_spin(-p, p-q+r-s, p-s)
	* occupation_numbers[p-q+r-s];
ns_value -= 8.000000 * current.interactions_same_spin(-q, -s, r)
	* sign(current.epsilon_tilde[-p] + current.epsilon_tilde[-r-s] - current.epsilon_tilde[-s] - current.epsilon_tilde[-p-r]) 
	* current.interactions_same_spin(-p, -r-s, p-s)
	* occupation_numbers[r+s];
ns_value -= 8.000000 * current.interactions_same_spin(p, q, r)
	* sign(current.epsilon_tilde[q] + current.epsilon_tilde[s] - current.epsilon_tilde[q] - current.epsilon_tilde[s]) 
	* current.interactions_differing_spin(q, s, Gamma<L>)
	* occupation_numbers[q];
dHdl.interactions_same_spin(p, q, r) += one_value + occupation_numbers[s] * ns_value;
} // s-loop
} // r-loop
} // q-loop
} // p-loop
dHdl.interactions_same_spin.symmetrize();
dHdl.interactions_differing_spin.symmetrize();
}
} // namespace NickelCUT::flow