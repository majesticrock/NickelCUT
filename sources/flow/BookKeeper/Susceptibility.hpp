#pragma once
#include "base.hpp"

#include "../../susceptibility/Charge.hpp"
#include "../../susceptibility/SC.hpp"
#include "../../susceptibility/Spin.hpp"

#include <nlohmann/json.hpp>

#include <chrono>
#include <deque>
#include <utility>
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

    std::deque<std::vector<double>> extremal_sc_eigenvectors;
    std::deque<std::vector<double>> extremal_charge_eigenvectors;
    std::deque<std::vector<double>> extremal_spin_eigenvectors;

    std::deque<std::pair<int, int>> extremal_sc_momenta;
    std::deque<std::pair<int, int>> extremal_charge_momenta;
    std::deque<std::pair<int, int>> extremal_spin_momenta;

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

    double initial_lambda_ex_derivative;
    // Guards against accidental 0-derivatives by requiring two adjacent derivatives to be nearly 0
    bool previous_was_stable = false;

    void compute_static_ieom();

    double compute_lambda_ex_derivative() noexcept;
};

}  // namespace flow::BookKeeper
}  // namespace NickelCUT
