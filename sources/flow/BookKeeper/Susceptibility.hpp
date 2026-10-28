#pragma once
#include "base.hpp"

#include "../../susceptibility/Charge.hpp"
#include "../../susceptibility/SC.hpp"
#include "../../susceptibility/Spin.hpp"

#include <nlohmann/json.hpp>

#include <chrono>
#include <deque>
#include <vector>

namespace NickelCUT {
struct Model;

namespace flow::BookKeeper {
struct Susceptibility : private base {
    typedef std::array<std::vector<double>, number_of_high_symmtry_momenta<L>> static_chi_container;

    std::deque<static_chi_container> sc_susceptibilities;
    std::deque<static_chi_container> charge_susceptibilities;
    std::deque<static_chi_container> spin_susceptibilities;

    std::deque<double> extremal_sc_eigenvalues;
    std::deque<double> extremal_charge_eigenvalues;
    std::deque<double> extremal_spin_eigenvalues;

    /////////////////////////////////////////////////////////

    Susceptibility(const FlowContainer& initial_flow_state,
                   FlowEquation& flow_equation,
                   double _band_width,
                   double _dl,
                   std::chrono::minutes::rep _max_runtime_duration,
                   const Model& _model);

    void print_final(const FlowContainer& x, double l);

    void operator()(const FlowContainer& x, double l);

    inline const FlowContainer& get_last_good_state() const noexcept {
        return last_good_state;
    }

    friend void to_json(nlohmann::json& j, const Susceptibility& book_keeper) noexcept;

private:
    const Model& model;

    susceptibility::SC ieom_sc;
    susceptibility::Charge ieom_charge;
    susceptibility::Spin ieom_spin;

    const std::vector<susceptibility::ieom_base::Vector> states;

    void compute_static_ieom();
};

}  // namespace flow::BookKeeper
}  // namespace NickelCUT
