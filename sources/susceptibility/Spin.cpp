#include "Spin.hpp"

#include "../flow/momentum_iterator.hpp"

namespace NickelCUT::susceptibility {
Spin::Spin(const flow::FlowContainer& _flow_result, const flow::Model& model) : ieom_base(_flow_result, model) {}

Spin::Spin(const ieom_base* const other_instance) : ieom_base(other_instance) {}

void Spin::fill_M(const momentum& x) {
    dynamical_matrix.setZero(N, N);
    
    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        const momentum kx = k + x;
        const double row_factor = difference_quotient(k, kx);

        dynamical_matrix(k.get_position(), k.get_position()) = -0.5 * row_factor;

        for (momentum l = momentum::begin(); l != momentum::end(); ++l) {
            const momentum lx = l + x;
            const double col_factor = difference_quotient(l, lx);
            dynamical_matrix(k.get_position(), l.get_position()) += (
                    - flow_result.interactions_differing_spin(k, lx, x) 
                    + flow_result.interactions_same_spin(k, lx, x)
                    - flow_result.interactions_same_spin(k, lx, l - k)
                ) * col_factor * row_factor;
        }
    }
}

void Spin::fill_N(const momentum& x) {
    norm_matrix.setZero(N);
    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        const momentum kx = k + x;
        norm_matrix(k.get_position()) = 0.5 * difference_quotient(k, kx);//2 * (occupation_numbers[k] - occupation_numbers[kx]);
    }
}

double Spin::difference_quotient(const momentum& k, const momentum& kx) const noexcept {
    const double energy_difference = dispersion[k] - dispersion[kx];
    if (is_zero(energy_difference)) {
        // Trivially, occupation_difference is now 0 as well, so we get the derivative of the Fermi function
        // at epsilon = 0, being -beta / 4;
        return -0.25 * beta;
    }
    return (occupation_numbers[k] - occupation_numbers[kx]) / energy_difference;
}

}  // namespace NickelCUT::susceptibility