#include "Spin.hpp"

#include "../momentum_iterator.hpp"

namespace NickelCUT::susceptibility {
Spin::Spin(const FlowContainer& _flow_result, const Model& model) : ieom_base(_flow_result, model) {}

Spin::Spin(const ieom_base* const other_instance) : ieom_base(other_instance) {}

void Spin::fill_U_matrix(const momentum& x) {
    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        for (momentum l = momentum::begin(); l != momentum::end(); ++l) {
            const momentum lx = l + x;
            U_matrix(k.get_position(), l.get_position()) = (
                    - flow_result->interactions_differing_spin(k, lx, x) 
                    + flow_result->interactions_same_spin(k, lx, x)
                    - flow_result->interactions_same_spin(k, lx, l - k)
                );
        }
    }
}

void Spin::fill_chi_matrix(const momentum& x) {
    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        chi_matrix(k.get_position()) = chi_number(k, k + x);

        if(chi_matrix(k.get_position()) < 0.) {
            std::cerr << "Spin: Negative chi_number!\n"
                << "Ek       = " << dispersion[k] << '\n'
                << "Ekx      = " << dispersion[k+x] << '\n'
                << "dE       = " << dispersion[k] - dispersion[k+x] << '\n'
                << "fk       = " << occupation_numbers[k] << '\n'
                << "fkx      = " << occupation_numbers[k+x] << '\n'
                << "df       = " << occupation_numbers[k+x] - occupation_numbers[k] << '\n'
                << "chi      = " << chi_matrix(k.get_position()) << '\n'
                << "beta     = " << beta << '\n';
        }
    }
}

}  // namespace NickelCUT::susceptibility