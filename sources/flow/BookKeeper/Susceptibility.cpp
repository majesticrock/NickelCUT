#include "Susceptibility.hpp"

#include "../ExtractionContainer.hpp"
#include "../../shared/Model.hpp"
#include "../../shared/FlowContainer.hpp"

#include <cmath>

namespace NickelCUT::flow::BookKeeper {
Susceptibility::Susceptibility(const FlowContainer& initial_flow_state,
                       FlowEquation& flow_equation,
                       double _band_width,
                       double _dl,
                       std::chrono::minutes::rep _max_runtime_duration,
                   const Model& _model)
    : base(initial_flow_state, flow_equation, _band_width, _dl, _max_runtime_duration),
    model(_model),
    ieom_sc(initial_flow_state, model),
    ieom_charge(&ieom_sc),
    ieom_spin(&ieom_sc),
    states({ieom_sc.swave_state(), ieom_sc.dwave_state()})
{
    compute_static_ieom();
}

void Susceptibility::print_final(const FlowContainer& x, double l) {
    if (basic_wrap_up(x, l)) {
        ieom_sc.update(x, model);
        ieom_charge.update(&ieom_sc);
        ieom_spin.update(&ieom_sc);
        compute_static_ieom();
    }
}

void Susceptibility::operator()(const FlowContainer& x, double l) {
    if (basic_processing(x, l)) {
        ieom_sc.update(x, model);
        ieom_charge.update(&ieom_sc);
        ieom_spin.update(&ieom_sc);
        compute_static_ieom();
    }
}

void Susceptibility::compute_static_ieom() {
    sc_susceptibilities.push_back(ieom_sc.static_susceptibility_at_high_symmetry(states));
    spin_susceptibilities.push_back(ieom_spin.static_susceptibility_at_high_symmetry(states));
    charge_susceptibilities.push_back(ieom_charge.static_susceptibility_at_high_symmetry(states));

    extremal_sc_eigenvalues.push_back(ieom_sc.get_extremal_eigenvalue());
    extremal_charge_eigenvalues.push_back(ieom_charge.get_extremal_eigenvalue());
    extremal_spin_eigenvalues.push_back(ieom_spin.get_extremal_eigenvalue());

    if (extremal_sc_eigenvalues.size() < 2U) return;

    const double sc_ev_shift = *(extremal_sc_eigenvalues.rbegin()) - *(extremal_sc_eigenvalues.rbegin() + 1);
    const double charge_ev_shift = *(extremal_charge_eigenvalues.rbegin()) - *(extremal_charge_eigenvalues.rbegin() + 1);
    const double spin_ev_shift = *(extremal_spin_eigenvalues.rbegin()) - *(extremal_spin_eigenvalues.rbegin() + 1);

    const double total_ev_shift = std::sqrt(
        sc_ev_shift * sc_ev_shift 
        + charge_ev_shift * charge_ev_shift 
        + spin_ev_shift * spin_ev_shift);
    
    std::cout << "Total eigenvalue shift is " << total_ev_shift << std::endl;

    // TODO throw exception for convergence or divergence.
}

void to_json(nlohmann::json& j, const Susceptibility& book_keeper) noexcept {
    j = book_keeper.convert_to_json();
    j.update({
        {"sc_susceptibilities", book_keeper.sc_susceptibilities},
        {"charge_susceptibilities", book_keeper.charge_susceptibilities},
        {"spin_susceptibilities", book_keeper.spin_susceptibilities},
        {"extremal_sc_eigenvalues", book_keeper.extremal_sc_eigenvalues},
        {"extremal_charge_eigenvalues", book_keeper.extremal_charge_eigenvalues},
        {"extremal_spin_eigenvalues", book_keeper.extremal_spin_eigenvalues},
        {"extracted_channels", ExtractionContainer(book_keeper.last_good_state)}
    });
}

}  // namespace NickelCUT::flow::BookKeeper
