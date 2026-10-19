#include "Spin.hpp"
#include "../flow/momentum_iterator.hpp"

namespace NickelCUT::susceptibility
{
Spin::Spin(const flow::FlowContainer& _flow_result, const flow::Model& model)
    : ieom_base(_flow_result, model)
{}

Spin::Spin(const ieom_base* const other_instance) 
    : ieom_base(other_instance) 
{}

void Spin::fill_M(const momentum& x) 
{
    dynamical_matrix.setZero(N, N);
    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        const momentum kx = k+x;
        double row_factor = 4. * (occupation_numbers[kx] - occupation_numbers[k]);

        dynamical_matrix(k.get_position(), k.get_position()) = 0.5 * (epsilon_tilde[k] - epsilon_tilde[kx]);
        for (momentum P = momentum::begin(); P != momentum::end(); ++P) {
            dynamical_matrix(k.get_position(), k.get_position()) -= (
              flow_result.interactions_differing_spin(kx, P, flow::Gamma<L>)
            - flow_result.interactions_differing_spin(k,  P, flow::Gamma<L>)
            - flow_result.interactions_same_spin(kx, P, P-kx)
            + flow_result.interactions_same_spin(kx, P, flow::Gamma<L>)
            + flow_result.interactions_same_spin(k,  P, P-k)
            - flow_result.interactions_same_spin(k,  P, flow::Gamma<L>)
        ) * occupation_numbers[P];
        }
        dynamical_matrix(k.get_position(), k.get_position()) *= row_factor;

        for (momentum l = momentum::begin(); l != momentum::end(); ++l) {
            const momentum lx = l+x;
            double col_factor = (occupation_numbers[lx] - occupation_numbers[l]);
            dynamical_matrix(k.get_position(), l.get_position()) += (
                - flow_result.interactions_differing_spin(k, lx, x)
                + flow_result.interactions_same_spin(     k, lx, x)
                - flow_result.interactions_same_spin(     k, lx, l-k)
            ) * col_factor * row_factor;
        }
    }
}

void Spin::fill_N(const momentum& x) 
{
    norm_matrix.setZero(N);
    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        const momentum kx = k+x;
        norm_matrix(k.get_position()) = 2 * (occupation_numbers[k] - occupation_numbers[kx]);
    }
}

} // namespace NickelCUT::susceptibility