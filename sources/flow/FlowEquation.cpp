#include "FlowEquation.hpp"

#include "momentum_iterator.hpp"
#include "FlowContainer.hpp"

namespace NickelCUT::flow {

void FlowEquation::operator()(const FlowContainer& current, FlowContainer& dHdl, const double /*l*/) {
for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
dHdl.dispersion[p.get_position()] = 0.0;
for (momentum_iterator<L> q = momentum_iterator<L>::begin(); q != momentum_iterator<L>::end(); ++q) {
double nr_value{};
double one_value{};
for (momentum_iterator<L> r = momentum_iterator<L>::begin(); r != momentum_iterator<L>::end(); ++r) {
one_value += current.interactions_same_spin((-r).get_position(), (-p+q+r).get_position(), (q+r).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(-r).get_position()] - current.epsilon_tilde[(-p+q+r).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (q).get_position(), (p-r).get_position());
one_value += current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position());
one_value += current.interactions_differing_spin((-r).get_position(), (-p+q+r).get_position(), (q+r).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(-p+q+r).get_position()] - current.epsilon_tilde[(-r).get_position()]) 
	* current.interactions_differing_spin((-p).get_position(), (q).get_position(), (q+r).get_position());
one_value += current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position());
one_value += current.interactions_differing_spin((-r).get_position(), (-p+q+r).get_position(), (q+r).get_position())
	* sign(current.epsilon_tilde[(-r).get_position()] + current.epsilon_tilde[(-p+q+r).get_position()] - current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(-p).get_position()]) 
	* current.interactions_differing_spin((-r).get_position(), (-p+q+r).get_position(), (q+r).get_position());
one_value += current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-q).get_position()] - current.epsilon_tilde[(-p).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (-q+r).get_position(), (p-q+r).get_position());
one_value += current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(-q).get_position()] - current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p-r).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (-q).get_position(), (p-q+r).get_position());
one_value += current.interactions_same_spin((-r).get_position(), (-p+q+r).get_position(), (q+r).get_position())
	* sign(current.epsilon_tilde[(r).get_position()] + current.epsilon_tilde[(p-q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((r).get_position(), (p-q-r).get_position(), (p-r).get_position());
nr_value += current.interactions_differing_spin((-p).get_position(), (p-q-r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(p-q-r).get_position()] - current.epsilon_tilde[(-q).get_position()] - current.epsilon_tilde[(-r).get_position()]) 
	* current.interactions_differing_spin((-p).get_position(), (p-q-r).get_position(), (p-q).get_position());
nr_value += current.interactions_same_spin((-q).get_position(), (-p+q+r).get_position(), (q+r).get_position())
	* sign(current.epsilon_tilde[(q).get_position()] + current.epsilon_tilde[(p-q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(-r).get_position()]) 
	* current.interactions_same_spin((q).get_position(), (p-q-r).get_position(), (p-q).get_position());
nr_value += current.interactions_differing_spin((p).get_position(), (q).get_position(), (q-r).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+q-r).get_position()] - current.epsilon_tilde[(r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (q-r).get_position());
nr_value += current.interactions_differing_spin((q).get_position(), (p-q-r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(r).get_position()] - current.epsilon_tilde[(-q).get_position()] - current.epsilon_tilde[(-p+q+r).get_position()]) 
	* current.interactions_differing_spin((-p).get_position(), (r).get_position(), (p-q).get_position());
nr_value += current.interactions_differing_spin((-q).get_position(), (-p+q+r).get_position(), (q+r).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(r).get_position()] - current.epsilon_tilde[(-p+q+r).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p).get_position(), (r).get_position(), (q+r).get_position());
nr_value += current.interactions_differing_spin((q).get_position(), (p-q-r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(q).get_position()] + current.epsilon_tilde[(p-q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(-r).get_position()]) 
	* current.interactions_differing_spin((q).get_position(), (p-q-r).get_position(), (p-q).get_position());
nr_value += current.interactions_differing_spin((-p).get_position(), (p-q-r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(q).get_position()] + current.epsilon_tilde[(r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(-p+q+r).get_position()]) 
	* current.interactions_differing_spin((q).get_position(), (r).get_position(), (p-q).get_position());
nr_value += current.interactions_differing_spin((-p).get_position(), (-r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(-r).get_position()] - current.epsilon_tilde[(-q).get_position()] - current.epsilon_tilde[(-p+q-r).get_position()]) 
	* current.interactions_differing_spin((-p).get_position(), (-r).get_position(), (p-q).get_position());
nr_value += current.interactions_differing_spin((q).get_position(), (r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(q).get_position()] + current.epsilon_tilde[(r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(-p+q+r).get_position()]) 
	* current.interactions_differing_spin((q).get_position(), (r).get_position(), (p-q).get_position());
nr_value += current.interactions_differing_spin((-p).get_position(), (-r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(q).get_position()] + current.epsilon_tilde[(p-q+r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(r).get_position()]) 
	* current.interactions_differing_spin((q).get_position(), (p-q+r).get_position(), (p-q).get_position());
nr_value += current.interactions_differing_spin((p).get_position(), (q).get_position(), (q-r).get_position())
	* sign(current.epsilon_tilde[(r).get_position()] + current.epsilon_tilde[(p+q-r).get_position()] - current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p).get_position()]) 
	* current.interactions_differing_spin((r).get_position(), (p+q-r).get_position(), (q-r).get_position());
nr_value += current.interactions_differing_spin((q).get_position(), (r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(p-q-r).get_position()] - current.epsilon_tilde[(-q).get_position()] - current.epsilon_tilde[(-r).get_position()]) 
	* current.interactions_differing_spin((-p).get_position(), (p-q-r).get_position(), (p-q).get_position());
nr_value += current.interactions_differing_spin((-q).get_position(), (-p+q+r).get_position(), (q+r).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(-p+q+r).get_position()] - current.epsilon_tilde[(r).get_position()] - current.epsilon_tilde[(-p).get_position()]) 
	* current.interactions_differing_spin((-q).get_position(), (-p+q+r).get_position(), (q+r).get_position());
nr_value += current.interactions_same_spin((q).get_position(), (r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(p-q-r).get_position()] - current.epsilon_tilde[(-r).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (p-q-r).get_position(), (p-r).get_position());
nr_value += current.interactions_same_spin((p).get_position(), (q).get_position(), (q-r).get_position())
	* sign(current.epsilon_tilde[(r).get_position()] + current.epsilon_tilde[(p+q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_same_spin((r).get_position(), (p+q-r).get_position(), (p-r).get_position());
nr_value += current.interactions_same_spin((-p).get_position(), (-q).get_position(), (p-r).get_position())
	* sign(current.epsilon_tilde[(r).get_position()] + current.epsilon_tilde[(p+q-r).get_position()] - current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p).get_position()]) 
	* current.interactions_same_spin((r).get_position(), (p+q-r).get_position(), (q-r).get_position());
nr_value += current.interactions_same_spin((p).get_position(), (q).get_position(), (q-r).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(-q).get_position()] - current.epsilon_tilde[(-r).get_position()] - current.epsilon_tilde[(-p-q+r).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (-q).get_position(), (p-r).get_position());
nr_value += current.interactions_same_spin((-p).get_position(), (-q).get_position(), (p-r).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+q-r).get_position()] - current.epsilon_tilde[(r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (q-r).get_position());
nr_value += current.interactions_same_spin((q).get_position(), (r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(-r).get_position()] - current.epsilon_tilde[(p-q-r).get_position()] - current.epsilon_tilde[(-p).get_position()]) 
	* current.interactions_same_spin((-q).get_position(), (-r).get_position(), (p-r).get_position());
nr_value += current.interactions_same_spin((-p).get_position(), (p-q-r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(-r).get_position()] - current.epsilon_tilde[(p-q-r).get_position()] - current.epsilon_tilde[(-p).get_position()]) 
	* current.interactions_same_spin((-q).get_position(), (-r).get_position(), (p-r).get_position());
nr_value += current.interactions_same_spin((q).get_position(), (p-q-r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(r).get_position()] - current.epsilon_tilde[(-p+q+r).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (r).get_position(), (q+r).get_position());
nr_value += current.interactions_same_spin((-q).get_position(), (-p+q+r).get_position(), (q+r).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(r).get_position()] - current.epsilon_tilde[(-q).get_position()] - current.epsilon_tilde[(-p+q+r).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (r).get_position(), (p-q).get_position());
nr_value += current.interactions_same_spin((-p).get_position(), (p-q-r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(p-q-r).get_position()] - current.epsilon_tilde[(-r).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (p-q-r).get_position(), (p-r).get_position());
nr_value += current.interactions_same_spin((q).get_position(), (p-q-r).get_position(), (p-q).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(-p+q+r).get_position()] - current.epsilon_tilde[(r).get_position()] - current.epsilon_tilde[(-p).get_position()]) 
	* current.interactions_same_spin((-q).get_position(), (-p+q+r).get_position(), (q+r).get_position());
nr_value *= current.occupation_numbers[r.get_position()];
} // r-loop
dHdl.dispersion[p.get_position()] += (nr_value + one_value) * current.occupation_numbers[q.get_position()];
} // q-loop
} // p-loop

//--------------------------------------------------------------//
for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
for (momentum_iterator<L> q = momentum_iterator<L>::begin(); q != momentum_iterator<L>::end(); ++q) {
for (momentum_iterator<L> r = momentum_iterator<L>::begin(); r != momentum_iterator<L>::end(); ++r) {
dHdl.interactions_differing_spin(p.get_position(), q.get_position(), r.get_position()) = 0.0;
dHdl.interactions_differing_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.epsilon_tilde[(q).get_position()];
dHdl.interactions_differing_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.epsilon_tilde[(q-r).get_position()];
dHdl.interactions_differing_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.epsilon_tilde[(p+r).get_position()];
dHdl.interactions_differing_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p+r).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_differing_spin((p+r).get_position(), (q-r).get_position(), (-r).get_position())
	* current.epsilon_tilde[(q-r).get_position()];
dHdl.interactions_differing_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p+r).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_differing_spin((p+r).get_position(), (q-r).get_position(), (-r).get_position())
	* current.epsilon_tilde[(p+r).get_position()];
dHdl.interactions_differing_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p+r).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_differing_spin((p+r).get_position(), (q-r).get_position(), (-r).get_position())
	* current.epsilon_tilde[(q).get_position()];
dHdl.interactions_differing_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.epsilon_tilde[(q).get_position()]
	* current.occupation_numbers[(q).get_position()];
dHdl.interactions_differing_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p+r).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_differing_spin((p+r).get_position(), (q-r).get_position(), (-r).get_position())
	* current.epsilon_tilde[(q).get_position()]
	* current.occupation_numbers[(q).get_position()];
dHdl.interactions_differing_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p+r).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_differing_spin((p+r).get_position(), (q-r).get_position(), (-r).get_position())
	* current.epsilon_tilde[(p).get_position()]
	* current.occupation_numbers[(p).get_position()];
dHdl.interactions_differing_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.epsilon_tilde[(p).get_position()]
	* current.occupation_numbers[(p).get_position()];
double ns_value{};
double one_value{};
for (momentum_iterator<L> s = momentum_iterator<L>::begin(); s != momentum_iterator<L>::end(); ++s) {
one_value += current.interactions_differing_spin((-s).get_position(), (p+q+s).get_position(), (p+r+s).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(-q).get_position()] - current.epsilon_tilde[(s).get_position()] - current.epsilon_tilde[(-p-q-s).get_position()]) 
	* current.interactions_differing_spin((-p).get_position(), (-q).get_position(), (p+s).get_position());
one_value += current.interactions_differing_spin((p).get_position(), (q).get_position(), (s).get_position())
	* sign(current.epsilon_tilde[(p+s).get_position()] + current.epsilon_tilde[(q-s).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p+s).get_position(), (q-s).get_position(), (r-s).get_position());
one_value += current.interactions_differing_spin((p).get_position(), (q).get_position(), (s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p-s).get_position()] - current.epsilon_tilde[(-q+s).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (-q+r).get_position(), (r-s).get_position());
one_value += current.interactions_differing_spin((-s).get_position(), (p+q+s).get_position(), (p+r+s).get_position())
	* sign(current.epsilon_tilde[(-s).get_position()] + current.epsilon_tilde[(p+q+s).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_differing_spin((-s).get_position(), (p+q+s).get_position(), (p+s).get_position());
ns_value += current.interactions_differing_spin((-s).get_position(), (p+q+s).get_position(), (p+r+s).get_position())
	* sign(current.epsilon_tilde[(-s).get_position()] + current.epsilon_tilde[(p+q+s).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_differing_spin((-s).get_position(), (p+q+s).get_position(), (p+s).get_position());
ns_value += current.interactions_same_spin((s).get_position(), (q-r).get_position(), (q-r-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((-p).get_position(), (p-q+r-s).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(p-q+r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (p-q+r-s).get_position(), (p+r-s).get_position());
ns_value += current.interactions_same_spin((s).get_position(), (q-r).get_position(), (q-r-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((p).get_position(), (r-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-q+r).get_position()] + current.epsilon_tilde[(-r+s).get_position()] - current.epsilon_tilde[(s).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-q+r).get_position(), (-r+s).get_position(), (q-r+s).get_position());
ns_value += current.interactions_differing_spin((-s).get_position(), (p+q+s).get_position(), (p+r+s).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(-q).get_position()] - current.epsilon_tilde[(s).get_position()] - current.epsilon_tilde[(-p-q-s).get_position()]) 
	* current.interactions_differing_spin((-p).get_position(), (-q).get_position(), (p+s).get_position());
ns_value += current.interactions_same_spin((s).get_position(), (p+r).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((-s).get_position(), (p+q+s).get_position(), (q-r+s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+q+s).get_position()] - current.epsilon_tilde[(-s).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (q+s).get_position());
ns_value += current.interactions_differing_spin((p).get_position(), (r-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(q-r).get_position()] + current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(-s).get_position()]) 
	* current.interactions_same_spin((q-r).get_position(), (r-s).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((p).get_position(), (r-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(s).get_position()] - current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-r+s).get_position()]) 
	* current.interactions_same_spin((-q).get_position(), (s).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((s).get_position(), (p+r).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p).get_position()]) 
	* current.interactions_differing_spin((s).get_position(), (p+r).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((s).get_position(), (q-r).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((q).get_position(), (s).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((s).get_position(), (q-r).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(s).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (r+s).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((q).get_position(), (s).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((s).get_position(), (p+r).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-p).get_position(), (-r+s).get_position(), (p+s).get_position())
	* sign(current.epsilon_tilde[(q-r).get_position()] + current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(-s).get_position()]) 
	* current.interactions_differing_spin((q-r).get_position(), (r-s).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-r-s).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p).get_position()]) 
	* current.interactions_differing_spin((s).get_position(), (p+r).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((-p).get_position(), (-q).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+q-s).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((s).get_position(), (p+q-s).get_position(), (p+r-s).get_position());
ns_value += current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(-p+q-r+s).get_position()]) 
	* current.interactions_differing_spin((s).get_position(), (q-r).get_position(), (p-s).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(r+s).get_position()]) 
	* current.interactions_same_spin((s).get_position(), (p+r).get_position(), (p-s).get_position());
ns_value += current.interactions_differing_spin((p).get_position(), (q).get_position(), (q-s).get_position())
	* sign(current.epsilon_tilde[(p+r).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p+q-s).get_position()] - current.epsilon_tilde[(s).get_position()]) 
	* current.interactions_differing_spin((p+r).get_position(), (q-r).get_position(), (q-r-s).get_position());
ns_value += current.interactions_differing_spin((-p).get_position(), (-q).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-p-q+s).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (-q+r).get_position(), (p+r-s).get_position());
ns_value += current.interactions_differing_spin((p).get_position(), (r-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(s).get_position()] - current.epsilon_tilde[(-r+s).get_position()] - current.epsilon_tilde[(-q+r).get_position()]) 
	* current.interactions_same_spin((-q).get_position(), (s).get_position(), (q-r+s).get_position());
ns_value += current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(-r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-p-r).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (-r-s).get_position(), (p-s).get_position());
ns_value += current.interactions_differing_spin((-p).get_position(), (p-q+r-s).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(p-q+r-s).get_position()] - current.epsilon_tilde[(-p-r).get_position()]) 
	* current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (p+r-s).get_position());
ns_value += current.interactions_same_spin((s).get_position(), (p+r).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((p).get_position(), (q).get_position(), (q-s).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+q-s).get_position()] - current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p+r).get_position()]) 
	* current.interactions_differing_spin((s).get_position(), (p+q-s).get_position(), (q-r-s).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(p-q+r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-q+r).get_position()]) 
	* current.interactions_differing_spin((-p).get_position(), (p-q+r-s).get_position(), (p-s).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-r-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(s).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (r+s).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-p).get_position(), (-r+s).get_position(), (p+s).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(s).get_position()] - current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-r+s).get_position()]) 
	* current.interactions_differing_spin((-q).get_position(), (s).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((-s).get_position(), (p+q+s).get_position(), (q-r+s).get_position())
	* sign(current.epsilon_tilde[(-s).get_position()] + current.epsilon_tilde[(p+q+s).get_position()] - current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p).get_position()]) 
	* current.interactions_differing_spin((-s).get_position(), (p+q+s).get_position(), (q+s).get_position());
ns_value += current.interactions_differing_spin((p).get_position(), (s).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position())
	* current.occupation_numbers[(p).get_position()];
ns_value += current.interactions_same_spin((q).get_position(), (s).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position())
	* current.occupation_numbers[(q).get_position()];
ns_value += current.interactions_differing_spin((p).get_position(), (s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-s).get_position()]) 
	* current.interactions_same_spin((-q).get_position(), (r-s).get_position(), (r).get_position())
	* current.occupation_numbers[(r-s).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(s).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (r+s).get_position(), (r).get_position())
	* current.occupation_numbers[(r+s).get_position()];
ns_value += current.interactions_same_spin((q).get_position(), (s).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.occupation_numbers[(q).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p).get_position()]) 
	* current.interactions_differing_spin((s).get_position(), (p+r).get_position(), (r).get_position())
	* current.occupation_numbers[(r+s).get_position()];
ns_value += current.interactions_differing_spin((p).get_position(), (s).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.occupation_numbers[(p).get_position()];
ns_value += current.interactions_differing_spin((p).get_position(), (s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-s).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-s).get_position(), (-q+r).get_position(), (r).get_position())
	* current.occupation_numbers[(r-s).get_position()];
ns_value += current.interactions_same_spin((-p).get_position(), (-s).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(-s).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-s).get_position(), (-q+r).get_position(), (r).get_position())
	* current.occupation_numbers[(r-s).get_position()];
ns_value += current.interactions_same_spin((-p).get_position(), (-s).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position())
	* current.occupation_numbers[(p).get_position()];
ns_value += current.interactions_differing_spin((-p).get_position(), (p-q+r-s).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(p-q+r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (p-q+r-s).get_position(), (p+r-s).get_position())
	* current.occupation_numbers[(-p+q-r+s).get_position()];
ns_value += current.interactions_differing_spin((p).get_position(), (s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(-r+s).get_position()]) 
	* current.interactions_same_spin((s).get_position(), (q-r).get_position(), (q-s).get_position())
	* current.occupation_numbers[(r-s).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position())
	* current.occupation_numbers[(q).get_position()];
ns_value += current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(r+s).get_position()]) 
	* current.interactions_same_spin((s).get_position(), (p+r).get_position(), (p-s).get_position())
	* current.occupation_numbers[(r+s).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-r-s).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p).get_position()]) 
	* current.interactions_differing_spin((s).get_position(), (p+r).get_position(), (r).get_position())
	* current.occupation_numbers[(r+s).get_position()];
ns_value += current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(-p+q-r+s).get_position()]) 
	* current.interactions_differing_spin((s).get_position(), (q-r).get_position(), (p-s).get_position())
	* current.occupation_numbers[(p-q+r-s).get_position()];
ns_value += current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(-r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-p-r).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (-r-s).get_position(), (p-s).get_position())
	* current.occupation_numbers[(r+s).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-r-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(s).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (r+s).get_position(), (r).get_position())
	* current.occupation_numbers[(r+s).get_position()];
ns_value += current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(p-q+r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-q+r).get_position()]) 
	* current.interactions_differing_spin((-p).get_position(), (p-q+r-s).get_position(), (p-s).get_position())
	* current.occupation_numbers[(p-q+r-s).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.occupation_numbers[(q).get_position()];
ns_value += current.interactions_same_spin((-p).get_position(), (-s).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-s).get_position()]) 
	* current.interactions_differing_spin((-q).get_position(), (r-s).get_position(), (r).get_position())
	* current.occupation_numbers[(r-s).get_position()];
ns_value += current.interactions_same_spin((-p).get_position(), (-s).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.occupation_numbers[(p).get_position()];
ns_value += current.interactions_differing_spin((-p).get_position(), (p-q+r-s).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(p-q+r-s).get_position()] - current.epsilon_tilde[(-p-r).get_position()]) 
	* current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (p+r-s).get_position())
	* current.occupation_numbers[(-p+q-r+s).get_position()];
ns_value += current.interactions_differing_spin((p).get_position(), (s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-q+r).get_position()]) 
	* current.interactions_same_spin((-q).get_position(), (r-s).get_position(), (q-s).get_position())
	* current.occupation_numbers[(r-s).get_position()];
dHdl.interactions_differing_spin(p.get_position(), q.get_position(), r.get_position()) += one_value + current.occupation_numbers[s.get_position()] * ns_value;
} // s-loop
} // r-loop
} // q-loop
} // p-loop

//--------------------------------------------------------------//
for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
for (momentum_iterator<L> q = momentum_iterator<L>::begin(); q != momentum_iterator<L>::end(); ++q) {
for (momentum_iterator<L> r = momentum_iterator<L>::begin(); r != momentum_iterator<L>::end(); ++r) {
dHdl.interactions_same_spin(p.get_position(), q.get_position(), r.get_position()) = 0.0;
dHdl.interactions_same_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.epsilon_tilde[(q).get_position()];
dHdl.interactions_same_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.epsilon_tilde[(q-r).get_position()];
dHdl.interactions_same_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.epsilon_tilde[(p+r).get_position()];
dHdl.interactions_same_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p+r).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_same_spin((p+r).get_position(), (q-r).get_position(), (-r).get_position())
	* current.epsilon_tilde[(q-r).get_position()];
dHdl.interactions_same_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p+r).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_same_spin((p+r).get_position(), (q-r).get_position(), (-r).get_position())
	* current.epsilon_tilde[(p+r).get_position()];
dHdl.interactions_same_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p+r).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_same_spin((p+r).get_position(), (q-r).get_position(), (-r).get_position())
	* current.epsilon_tilde[(q).get_position()];
dHdl.interactions_same_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.epsilon_tilde[(q).get_position()]
	* current.occupation_numbers[(q).get_position()];
dHdl.interactions_same_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p+r).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_same_spin((p+r).get_position(), (q-r).get_position(), (-r).get_position())
	* current.epsilon_tilde[(q).get_position()]
	* current.occupation_numbers[(q).get_position()];
dHdl.interactions_same_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p+r).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_same_spin((p+r).get_position(), (q-r).get_position(), (-r).get_position())
	* current.epsilon_tilde[(p).get_position()]
	* current.occupation_numbers[(p).get_position()];
dHdl.interactions_same_spin(p.get_position(), q.get_position(), r.get_position()) += sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.epsilon_tilde[(p).get_position()]
	* current.occupation_numbers[(p).get_position()];
double ns_value{};
double one_value{};
for (momentum_iterator<L> s = momentum_iterator<L>::begin(); s != momentum_iterator<L>::end(); ++s) {
one_value += current.interactions_same_spin((-s).get_position(), (p+q+s).get_position(), (p+r+s).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(-q).get_position()] - current.epsilon_tilde[(s).get_position()] - current.epsilon_tilde[(-p-q-s).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (-q).get_position(), (p+s).get_position());
one_value += current.interactions_same_spin((p).get_position(), (q).get_position(), (s).get_position())
	* sign(current.epsilon_tilde[(p+s).get_position()] + current.epsilon_tilde[(q-s).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p+s).get_position(), (q-s).get_position(), (r-s).get_position());
one_value += current.interactions_same_spin((p).get_position(), (q).get_position(), (s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p-s).get_position()] - current.epsilon_tilde[(-q+s).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (-q+r).get_position(), (r-s).get_position());
one_value += current.interactions_same_spin((-s).get_position(), (p+q+s).get_position(), (p+r+s).get_position())
	* sign(current.epsilon_tilde[(-s).get_position()] + current.epsilon_tilde[(p+q+s).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_same_spin((-s).get_position(), (p+q+s).get_position(), (p+s).get_position());
ns_value += current.interactions_same_spin((-s).get_position(), (p+q+s).get_position(), (p+r+s).get_position())
	* sign(current.epsilon_tilde[(-s).get_position()] + current.epsilon_tilde[(p+q+s).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(q).get_position()]) 
	* current.interactions_same_spin((-s).get_position(), (p+q+s).get_position(), (p+s).get_position());
ns_value += current.interactions_same_spin((s).get_position(), (q-r).get_position(), (q-r-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-p).get_position(), (p-q+r-s).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(p-q+r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (p-q+r-s).get_position(), (p+r-s).get_position());
ns_value += current.interactions_same_spin((s).get_position(), (q-r).get_position(), (q-r-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((p).get_position(), (r-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-q+r).get_position()] + current.epsilon_tilde[(-r+s).get_position()] - current.epsilon_tilde[(s).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-q+r).get_position(), (-r+s).get_position(), (q-r+s).get_position());
ns_value += current.interactions_same_spin((-s).get_position(), (p+q+s).get_position(), (p+r+s).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(-q).get_position()] - current.epsilon_tilde[(s).get_position()] - current.epsilon_tilde[(-p-q-s).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (-q).get_position(), (p+s).get_position());
ns_value += current.interactions_same_spin((s).get_position(), (p+r).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-s).get_position(), (p+q+s).get_position(), (q-r+s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+q+s).get_position()] - current.epsilon_tilde[(-s).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (q+s).get_position());
ns_value += current.interactions_differing_spin((p).get_position(), (r-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(q-r).get_position()] + current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(-s).get_position()]) 
	* current.interactions_differing_spin((q-r).get_position(), (r-s).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((p).get_position(), (r-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(s).get_position()] - current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-r+s).get_position()]) 
	* current.interactions_differing_spin((-q).get_position(), (s).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((s).get_position(), (p+r).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p).get_position()]) 
	* current.interactions_differing_spin((s).get_position(), (p+r).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((s).get_position(), (q-r).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((q).get_position(), (s).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((s).get_position(), (q-r).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(s).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (r+s).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((q).get_position(), (s).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position());
ns_value += current.interactions_differing_spin((s).get_position(), (p+r).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-p).get_position(), (-r+s).get_position(), (p+s).get_position())
	* sign(current.epsilon_tilde[(q-r).get_position()] + current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(-s).get_position()]) 
	* current.interactions_same_spin((q-r).get_position(), (r-s).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-r-s).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p).get_position()]) 
	* current.interactions_same_spin((s).get_position(), (p+r).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-p).get_position(), (-q).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+q-s).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((s).get_position(), (p+q-s).get_position(), (p+r-s).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(-p+q-r+s).get_position()]) 
	* current.interactions_same_spin((s).get_position(), (q-r).get_position(), (p-s).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(r+s).get_position()]) 
	* current.interactions_same_spin((s).get_position(), (p+r).get_position(), (p-s).get_position());
ns_value += current.interactions_same_spin((p).get_position(), (q).get_position(), (q-s).get_position())
	* sign(current.epsilon_tilde[(p+r).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p+q-s).get_position()] - current.epsilon_tilde[(s).get_position()]) 
	* current.interactions_same_spin((p+r).get_position(), (q-r).get_position(), (q-r-s).get_position());
ns_value += current.interactions_same_spin((-p).get_position(), (-q).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-p-q+s).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (-q+r).get_position(), (p+r-s).get_position());
ns_value += current.interactions_same_spin((p).get_position(), (r-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(s).get_position()] - current.epsilon_tilde[(-r+s).get_position()] - current.epsilon_tilde[(-q+r).get_position()]) 
	* current.interactions_same_spin((-q).get_position(), (s).get_position(), (q-r+s).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(-r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-p-r).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (-r-s).get_position(), (p-s).get_position());
ns_value += current.interactions_same_spin((-p).get_position(), (p-q+r-s).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(p-q+r-s).get_position()] - current.epsilon_tilde[(-p-r).get_position()]) 
	* current.interactions_same_spin((-q).get_position(), (-s).get_position(), (p+r-s).get_position());
ns_value += current.interactions_same_spin((s).get_position(), (p+r).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((p).get_position(), (q).get_position(), (q-s).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+q-s).get_position()] - current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p+r).get_position()]) 
	* current.interactions_same_spin((s).get_position(), (p+q-s).get_position(), (q-r-s).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(p-q+r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-q+r).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (p-q+r-s).get_position(), (p-s).get_position());
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-r-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(s).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (r+s).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-p).get_position(), (-r+s).get_position(), (p+s).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(s).get_position()] - current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-r+s).get_position()]) 
	* current.interactions_same_spin((-q).get_position(), (s).get_position(), (r).get_position());
ns_value += current.interactions_same_spin((-s).get_position(), (p+q+s).get_position(), (q-r+s).get_position())
	* sign(current.epsilon_tilde[(-s).get_position()] + current.epsilon_tilde[(p+q+s).get_position()] - current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p).get_position()]) 
	* current.interactions_same_spin((-s).get_position(), (p+q+s).get_position(), (q+s).get_position());
ns_value += current.interactions_differing_spin((p).get_position(), (s).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position())
	* current.occupation_numbers[(p).get_position()];
ns_value += current.interactions_differing_spin((q).get_position(), (s).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position())
	* current.occupation_numbers[(q).get_position()];
ns_value += current.interactions_differing_spin((p).get_position(), (s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-s).get_position()]) 
	* current.interactions_differing_spin((-q).get_position(), (r-s).get_position(), (r).get_position())
	* current.occupation_numbers[(r-s).get_position()];
ns_value += current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(s).get_position()]) 
	* current.interactions_differing_spin((p).get_position(), (r+s).get_position(), (r).get_position())
	* current.occupation_numbers[(r+s).get_position()];
ns_value += current.interactions_differing_spin((q).get_position(), (s).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.occupation_numbers[(q).get_position()];
ns_value += current.interactions_differing_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p).get_position()]) 
	* current.interactions_differing_spin((s).get_position(), (p+r).get_position(), (r).get_position())
	* current.occupation_numbers[(r+s).get_position()];
ns_value += current.interactions_differing_spin((p).get_position(), (s).get_position(), GammaPoint<L>.get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.occupation_numbers[(p).get_position()];
ns_value += current.interactions_differing_spin((p).get_position(), (s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-s).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_differing_spin((-s).get_position(), (-q+r).get_position(), (r).get_position())
	* current.occupation_numbers[(r-s).get_position()];
ns_value += current.interactions_same_spin((-p).get_position(), (-s).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(-s).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-s).get_position(), (-q+r).get_position(), (r).get_position())
	* current.occupation_numbers[(r-s).get_position()];
ns_value += current.interactions_same_spin((-p).get_position(), (-s).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position())
	* current.occupation_numbers[(p).get_position()];
ns_value += current.interactions_same_spin((-p).get_position(), (p-q+r-s).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(p-q+r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (p-q+r-s).get_position(), (p+r-s).get_position())
	* current.occupation_numbers[(-p+q-r+s).get_position()];
ns_value += current.interactions_same_spin((p).get_position(), (s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(-r+s).get_position()]) 
	* current.interactions_same_spin((s).get_position(), (q-r).get_position(), (q-s).get_position())
	* current.occupation_numbers[(r-s).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-s).get_position())
	* sign(current.epsilon_tilde[(-p-r).get_position()] + current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-p).get_position()] - current.epsilon_tilde[(-q).get_position()]) 
	* current.interactions_same_spin((-p-r).get_position(), (-q+r).get_position(), (r).get_position())
	* current.occupation_numbers[(q).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(r+s).get_position()]) 
	* current.interactions_same_spin((s).get_position(), (p+r).get_position(), (p-s).get_position())
	* current.occupation_numbers[(r+s).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-r-s).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p).get_position()]) 
	* current.interactions_same_spin((s).get_position(), (p+r).get_position(), (r).get_position())
	* current.occupation_numbers[(r+s).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(s).get_position()] + current.epsilon_tilde[(q-r).get_position()] - current.epsilon_tilde[(p).get_position()] - current.epsilon_tilde[(-p+q-r+s).get_position()]) 
	* current.interactions_same_spin((s).get_position(), (q-r).get_position(), (p-s).get_position())
	* current.occupation_numbers[(p-q+r-s).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(-r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-p-r).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (-r-s).get_position(), (p-s).get_position())
	* current.occupation_numbers[(r+s).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-r-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(r+s).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(s).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (r+s).get_position(), (r).get_position())
	* current.occupation_numbers[(r+s).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(-p).get_position()] + current.epsilon_tilde[(p-q+r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-q+r).get_position()]) 
	* current.interactions_same_spin((-p).get_position(), (p-q+r-s).get_position(), (p-s).get_position())
	* current.occupation_numbers[(p-q+r-s).get_position()];
ns_value += current.interactions_same_spin((-q).get_position(), (-s).get_position(), (q-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.occupation_numbers[(q).get_position()];
ns_value += current.interactions_same_spin((-p).get_position(), (-s).get_position(), (p+r-s).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(-q+r).get_position()] - current.epsilon_tilde[(-s).get_position()]) 
	* current.interactions_same_spin((-q).get_position(), (r-s).get_position(), (r).get_position())
	* current.occupation_numbers[(r-s).get_position()];
ns_value += current.interactions_same_spin((-p).get_position(), (-s).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(p).get_position()] + current.epsilon_tilde[(q).get_position()] - current.epsilon_tilde[(p+r).get_position()] - current.epsilon_tilde[(q-r).get_position()]) 
	* current.interactions_same_spin((p).get_position(), (q).get_position(), (r).get_position())
	* current.occupation_numbers[(p).get_position()];
ns_value += current.interactions_same_spin((-p).get_position(), (p-q+r-s).get_position(), (p-s).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(p-q+r-s).get_position()] - current.epsilon_tilde[(-p-r).get_position()]) 
	* current.interactions_same_spin((-q).get_position(), (-s).get_position(), (p+r-s).get_position())
	* current.occupation_numbers[(-p+q-r+s).get_position()];
ns_value += current.interactions_same_spin((p).get_position(), (s).get_position(), (r).get_position())
	* sign(current.epsilon_tilde[(-q).get_position()] + current.epsilon_tilde[(r-s).get_position()] - current.epsilon_tilde[(-s).get_position()] - current.epsilon_tilde[(-q+r).get_position()]) 
	* current.interactions_same_spin((-q).get_position(), (r-s).get_position(), (q-s).get_position())
	* current.occupation_numbers[(r-s).get_position()];
dHdl.interactions_same_spin(p.get_position(), q.get_position(), r.get_position()) += one_value + current.occupation_numbers[s.get_position()] * ns_value;
} // s-loop
} // r-loop
} // q-loop
} // p-loop
}
} // namespace NickelCUT::flow