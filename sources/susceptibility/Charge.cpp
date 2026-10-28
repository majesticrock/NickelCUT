#include "Charge.hpp"

#include "../momentum_iterator.hpp"

namespace NickelCUT::susceptibility {
Charge::Charge(const FlowContainer& _flow_result, const Model& model) : ieom_base(_flow_result, model) {}

Charge::Charge(const ieom_base* const other_instance) : ieom_base(other_instance) {}

void Charge::fill_U_matrix(const momentum& x) {
    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        for (momentum l = momentum::begin(); l != momentum::end(); ++l) {
            const momentum lx = l + x;
            U_matrix(k.get_position(), l.get_position()) = (
                    + flow_result->interactions_differing_spin(k, lx, x) 
                    + flow_result->interactions_same_spin(k, lx, x)
                    - flow_result->interactions_same_spin(k, lx, l - k)
                );
        }
    }
}

void Charge::fill_chi_matrix(const momentum& x) {
    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        chi_matrix(k.get_position()) = chi_number(k, k + x);
    }
}

}  // namespace NickelCUT::susceptibility