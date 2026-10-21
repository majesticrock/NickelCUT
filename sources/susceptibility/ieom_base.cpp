#include "ieom_base.hpp"

#include "../L.hpp"
#include "../flow/DecouplingChannel.hpp"

#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>

#include <functional>
#include <utility>

namespace NickelCUT::susceptibility {
template <typename F>
struct selfconsistency_solver_interface {
    F compute;

    std::string info() const { return "Hartree-Fock for iEoM"; }

    void iteration_step(const ieom_base::Vector& init, ieom_base::Vector& result) { compute(init, result); }
};

//--------------------------------------------------------------------------------//

ieom_base::ieom_base(const flow::FlowContainer& _flow_result, const flow::Model& model)
    : dynamical_matrix{Matrix::Zero(N, N)},
      norm_matrix{Matrix::Zero(N, N)},
      flow_result{_flow_result},
      epsilon_tilde{flow_result.epsilon_tilde()},
      beta{model.beta} {
    solve_hartree_fock(model);
}

ieom_base::ieom_base(const ieom_base* const other_instance)
    : dynamical_matrix{Matrix::Zero(N, N)},
      norm_matrix{Matrix::Zero(N, N)},
      occupation_numbers{other_instance->occupation_numbers},
      dispersion{other_instance->dispersion},
      flow_result{other_instance->flow_result},
      epsilon_tilde{other_instance->epsilon_tilde},
      beta{other_instance->beta} {};

void ieom_base::fill_matrices(const momentum& x) {
    fill_M(x);
    fill_N(x);
}

std::vector<double> ieom_base::compute_static_susceptibility(const momentum& x, const std::vector<Vector>& states) {
    prepare_static(x);

    std::vector<double> result;
    result.reserve(states.size());
    for (const auto& state : states) {
        const Vector tmp = eigen_solver.eigenvectors().adjoint() * norm_matrix.asDiagonal() * state;
        Vector inverse_eigenvalues = eigen_solver.eigenvalues();
        for (auto& ev : inverse_eigenvalues) {
            if (is_finite(ev)) {
                ev = 1. / ev;
                if (ev < 0.0) {
                    _found_negative_excitation = true;
                }
            } else {
                ev = 0.;
            }
        }

        result.push_back(tmp.dot(inverse_eigenvalues.asDiagonal() * tmp));
    }

    return result;
}

std::array<std::vector<double>, flow::number_of_high_symmtry_momenta<L>>
ieom_base::static_susceptibility_at_high_symmetry(const std::vector<Vector>& states) {
    const auto momenta = flow::square_lattice_high_symmetry_path<L>();
    std::array<std::vector<double>, flow::number_of_high_symmtry_momenta<L>> result;

    for (std::size_t i = 0; i < flow::number_of_high_symmtry_momenta<L>; ++i) {
        result[i] = compute_static_susceptibility(momenta[i], states);
    }

    return result;
}

ieom_base::Vector ieom_base::swave_state() const {
    return Vector::Ones(N).normalized();
}

ieom_base::Vector ieom_base::dwave_state() const {
    Vector result = Vector::Zero(N);
    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        result(k.get_position()) = std::cos(k.get_kx()) - std::cos(k.get_ky());
    }
    result.normalize();
    return result;
}

const std::array<double, N>& ieom_base::get_dispersion() const noexcept {
    return dispersion;
}

bool ieom_base::found_negative_excitation() const noexcept {
    return _found_negative_excitation;
}

//-----------------------------------------------------//
// Private methods

void ieom_base::solve_hartree_fock(const flow::Model& model) {
    const auto compute_occupation_numbers = [this, &model](const Vector& dispersion) {
        for (int i = 0; i < N; ++i) {
            occupation_numbers[i] = model.fermi_function(dispersion[i]);
        }
    };

    Eigen::Map<const Vector> epsilon_tilde_map(epsilon_tilde.data(), N);
    Vector self_energy = Eigen::Map<const Vector>(flow_result.dispersion.data(), N) - epsilon_tilde_map;

    const auto [differing, same] = flow::DecouplingChannel::SingleParticleEnergy(flow_result);

    const auto compute_self_energy = [&](const Vector& init, Vector& result) {
        result.setZero(N);
        compute_occupation_numbers(epsilon_tilde_map + init);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                // 1/N is contained in the interaction itself, and we assume spin invariance
                result(i) += (differing(i, j) + same(i, j)) * occupation_numbers[j];
            }
        }
        self_energy = 0.5 * (self_energy + result);
        result -= self_energy;
    };
    selfconsistency_solver_interface interface(compute_self_energy);

    auto selfconsistency_solver =
        mrock::utility::Selfconsistency::make_broyden<double>(&interface, &self_energy, 400, 1e-6);
    mrock::utility::Selfconsistency::ConvergenceInfo<double> convergence =
        selfconsistency_solver.compute(true, 600, 1e-10);

    if (!convergence.converged) {
        std::cout << "Hartree-Fock calculations did not converge. Retrying..." << std::endl;
        convergence = selfconsistency_solver.compute(true, 600, 1e-10);

        if (!convergence.converged) {
            throw std::runtime_error("Hartree-Fock calculations did not converge even after retry.");
        }
    }
    compute_occupation_numbers(epsilon_tilde_map + self_energy);
    for (int i = 0; i < N; ++i) {
        dispersion[i] = epsilon_tilde[i] + self_energy(i);
    }
}

void ieom_base::prepare_static(const momentum& x) {
    fill_matrices(x);
    eigen_solver.compute(dynamical_matrix);
}

ChiContainer::ChiContainer(const Eigen::VectorXd& eigen_energies, std::size_t n_states)
    : energies(eigen_energies.data(), eigen_energies.data() + eigen_energies.size()) {
    weights.reserve(n_states);
}

void to_json(nlohmann::json& j, const ChiContainer& container) noexcept {
    j = nlohmann::json{{"energies", container.energies}, {"weights", container.weights}};
}

}  // namespace NickelCUT::susceptibility