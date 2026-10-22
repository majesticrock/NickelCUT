#pragma once

#include "../L.hpp"
#include "../flow/FlowContainer.hpp"
#include "../flow/Model.hpp"
#include "../flow/momentum_iterator.hpp"

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
    using momentum = flow::momentum_iterator<L>;

    ieom_base(const flow::FlowContainer& _flow_result, const flow::Model& model);
    ieom_base(const ieom_base* const other_instance);

    virtual ~ieom_base() = default;

    virtual void fill_chi_matrix(const momentum& x) = 0;
    virtual void fill_U_matrix(const momentum& x) = 0;

    void fill_matrices(const momentum& x);

    // Computes the static susceptibility chi(x;0) = <state | N [M]^{-1} N |state> for each state in states
    std::vector<double> compute_static_susceptibility(const momentum& x, const std::vector<Vector>& states);

    std::array<std::vector<double>, flow::number_of_high_symmtry_momenta<L>> static_susceptibility_at_high_symmetry(
        const std::vector<Vector>& states);

    Vector swave_state() const;
    Vector dwave_state() const;

    const std::array<double, N>& get_dispersion() const noexcept;
    bool found_negative_excitation() const noexcept;

protected:
    Matrix U_matrix;
    Vector chi_matrix; // This is diagonal in k, and we represent it as a vector
    Eigen::SelfAdjointEigenSolver<Matrix> eigen_solver;

    std::array<double, N> occupation_numbers;
    std::array<double, N> dispersion;

    const flow::FlowContainer& flow_result;
    const std::array<double, N> epsilon_tilde;

    const double beta;

    bool _found_negative_excitation = false;

    // This is positive definite!
    double chi_number(const momentum& k, const momentum& kx) const noexcept;
    // This is positive definite!
    double chi_pair(const momentum&k, const momentum& kx) const noexcept;

private:
    void solve_hartree_fock(const flow::Model& model);
    // Computes the Kernel K and diagoanlizes it, storing the result in eigen_solver
    void prepare_static(const momentum& x);
};

}  // namespace NickelCUT::susceptibility
