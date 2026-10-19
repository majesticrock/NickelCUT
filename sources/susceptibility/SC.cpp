#include "SC.hpp"
#include "../flow/momentum_iterator.hpp"

namespace NickelCUT::susceptibility
{
SC::SC(const flow::FlowContainer& _flow_result, const flow::Model& model)
    : ieom_base(_flow_result, model)
{}

SC::SC(const ieom_base* const other_instance) 
    : ieom_base(other_instance) 
{}

void SC::fill_M(const momentum& x) 
{
    dynamical_matrix.setZero(N, N);
    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        const momentum kx = k+x;
        double row_factor = 2. * (1. - occupation_numbers[k] - occupation_numbers[kx]);

        dynamical_matrix(k.get_position(), k.get_position()) = 0.5 * (epsilon_tilde[k] + epsilon_tilde[kx]);
        for (momentum P = momentum::begin(); P != momentum::end(); ++P) {
            dynamical_matrix(k.get_position(), k.get_position()) += (
                  flow_result.interactions_differing_spin(kx, P, flow::Gamma<L>)
                + flow_result.interactions_differing_spin(k,  P, flow::Gamma<L>)
                - flow_result.interactions_same_spin(kx, P, P-kx)
                + flow_result.interactions_same_spin(kx, P, flow::Gamma<L>)
                - flow_result.interactions_same_spin(k,  P, P-k)
                + flow_result.interactions_same_spin(k,  P, flow::Gamma<L>)
            ) * occupation_numbers[P];
        }
        dynamical_matrix(k.get_position(), k.get_position()) *= row_factor;

        for (momentum l = momentum::begin(); l != momentum::end(); ++l) {
            double col_factor = (1. - occupation_numbers[l] - occupation_numbers[l+x]);
            dynamical_matrix(k.get_position(), l.get_position()) += flow_result.interactions_differing_spin(-kx,k,k-l) 
                * col_factor * row_factor;
        }
    }
}

void SC::fill_N(const momentum& x) 
{
    norm_matrix.setZero(N);
    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        const momentum kx = k+x;
        norm_matrix(k.get_position()) = -(1. - occupation_numbers[k] - occupation_numbers[kx]);
    }
}

} // namespace NickelCUT::susceptibility