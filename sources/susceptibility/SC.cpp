#include "SC.hpp"

#include "../flow/momentum_iterator.hpp"

namespace NickelCUT::susceptibility {
SC::SC(const flow::FlowContainer& _flow_result, const flow::Model& model) : ieom_base(_flow_result, model) {}

SC::SC(const ieom_base* const other_instance) : ieom_base(other_instance) {}

void SC::fill_U_matrix(const momentum& x) {
    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        for (momentum l = momentum::begin(); l != momentum::end(); ++l) {
            const momentum lx = l + x;
            U_matrix(k.get_position(), l.get_position()) = 2. * flow_result.interactions_differing_spin(k, lx, x);
        }
    }
}

void SC::fill_chi_matrix(const momentum& x) {
    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        chi_matrix(k.get_position()) = chi_pair(k, k + x);

        if(chi_matrix(k.get_position()) < 0.) {
            std::cerr << "SC: Negative chi_number!\n"
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