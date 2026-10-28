#pragma once

#include "../L.hpp"
#include "../shared/FlowContainer.hpp"
#include "../shared/Model.hpp"
#include "../momentum_iterator.hpp"

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#include <array>
#include <vector>

namespace NickelCUT::susceptibility {

struct ChiContainer {
    std::vector<double> energies;
    std::vector<std::vector<double>> weights;

    ChiContainer() = default;
    ChiContainer(const Eigen::VectorXd& eigen_energies, std::size_t n_states);
};

void to_json(nlohmann::json& j, const ChiContainer& container) noexcept;

class ieom_base {
public:
    using Matrix = Eigen::MatrixXd;
    using Vector = Eigen::VectorXd;
    using momentum = momentum_iterator<L>;

    ieom_base(const FlowContainer& _flow_result, const Model& model);
    ieom_base(const ieom_base* const other_instance);

    void update(const FlowContainer& _flow_result, const Model& model);
    void update(const ieom_base* const other_instance);

    virtual ~ieom_base() = default;

    virtual void fill_chi_matrix(const momentum& x) = 0;
    virtual void fill_U_matrix(const momentum& x) = 0;

    void fill_matrices(const momentum& x);

    // Computes the static susceptibility chi(x;0) = <state | N [M]^{-1} N |state> for each state in states
    std::vector<double> compute_static_susceptibility(const momentum& x, const std::vector<Vector>& states);

    std::array<std::vector<double>, number_of_high_symmtry_momenta<L>> static_susceptibility_at_high_symmetry(
        const std::vector<Vector>& states);

    Vector swave_state() const;
    Vector dwave_state() const;

    const std::array<double, N>& get_dispersion() const noexcept;
    bool found_negative_excitation() const noexcept;
    bool found_pomeranchuk_instability() const noexcept;
    double get_extremal_eigenvalue() const noexcept;

    double hartree_fock_filling() const noexcept;

protected:
    Matrix U_matrix;
    Vector chi_matrix; // This is diagonal in k, and we represent it as a vector
    Eigen::SelfAdjointEigenSolver<Matrix> eigen_solver;

    std::array<double, N> occupation_numbers;
    std::array<double, N> dispersion;

    FlowContainer const * flow_result;
    std::array<double, N> epsilon_tilde;

    const double beta;
    double chemical_potential = 0.0;

    double _extremal_eigenvalue = 10000.;
    bool _found_negative_excitation = false;
    bool _found_pomeranchuk = false;

    // This is positive definite!
    double chi_number(const momentum& k, const momentum& kx) const noexcept;
    // This is positive definite!
    double chi_pair(const momentum&k, const momentum& kx) const noexcept;

private:
    double compute_chemical_potential(const Model& model, const double target_filling);
    void solve_hartree_fock(const Model& model);
    // Computes the Kernel K and diagoanlizes it, storing the result in eigen_solver
    void prepare_static(const momentum& x);
};

}  // namespace NickelCUT::susceptibility
