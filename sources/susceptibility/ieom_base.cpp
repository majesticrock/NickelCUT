#include "ieom_base.hpp"

#include "../shared/L.hpp"
#include "../flow/DecouplingChannel.hpp"

#include <boost/math/tools/toms748_solve.hpp>
#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <utility>

namespace NickelCUT::susceptibility {
template <typename F>
struct selfconsistency_solver_interface {
    F compute;

    std::string info() const { return "Hartree-Fock for iEoM"; }

    void iteration_step(const ieom_base::Vector& init, ieom_base::Vector& result) { compute(init, result); }
};

constexpr double log2 = 0.693147180559945309417232121458;

// Stable log(cosh(x)) for real x.
double log_cosh(double x) noexcept
{
    const double ax = std::abs(x);
    return ax + std::log1p(std::exp(-2.0 * ax)) - log2;
}

// Stable log(sinh(x) / x), for real x.
double log_sinhc(double x) noexcept
{
    const double ax = std::abs(x);
    if (ax < 1e-3) {
        // log(sinh(x)/x) = x²/6 - x⁴/180 + x⁶/2835 + O(x⁸)
        const double x2 = x * x;
        return x2 * (1.0 / 6.0
               + x2 * (-1.0 / 180.0
               + x2 / 2835.0));
    }
    // log(sinh(ax)/ax), written to avoid sinh overflow.
    return ax + std::log1p(-std::exp(-2.0 * ax))
               - std::log(2.0 * ax);
}

//--------------------------------------------------------------------------------//

ieom_base::ieom_base(const flow::FlowContainer& _flow_result, const Model& model)
    : U_matrix{Matrix::Zero(N, N)},
      chi_matrix{Vector::Zero(N)},
      flow_result{_flow_result},
      epsilon_tilde{flow_result.epsilon_tilde()},
      beta{model.beta} {
    solve_hartree_fock(model);
}

ieom_base::ieom_base(const ieom_base* const other_instance)
    : U_matrix{Matrix::Zero(N, N)},
      chi_matrix{Vector::Zero(N)},
      occupation_numbers{other_instance->occupation_numbers},
      dispersion{other_instance->dispersion},
      flow_result{other_instance->flow_result},
      epsilon_tilde{other_instance->epsilon_tilde},
      beta{other_instance->beta} {};

void ieom_base::fill_matrices(const momentum& x) {
    fill_U_matrix(x);
    fill_chi_matrix(x);
}

std::vector<double> ieom_base::compute_static_susceptibility(const momentum& x, const std::vector<Vector>& states) {
    prepare_static(x);

    std::vector<double> result;
    result.reserve(states.size());
    for (const auto& state : states) {
        const Vector tmp = eigen_solver.eigenvectors().adjoint() * chi_matrix.asDiagonal() * state;
        Vector inverse_eigenvalues = eigen_solver.eigenvalues();
        for (auto& ev : inverse_eigenvalues) {
            ev += 1.; // The matrix to invert is 1+K
            if (is_non_zero(ev)) {
                ev = 1. / ev;
                if (ev < 0.0) {
                    _found_negative_excitation = true;
                }
            } 
            else {
                std::cerr << "Warning: Inversion fallback was used!" << std::endl;
                ev = PRECISION;
            }
        }

        result.push_back(tmp.dot(inverse_eigenvalues.asDiagonal() * tmp));
    }

    return result;
}

std::array<std::vector<double>, number_of_high_symmtry_momenta<L>>
ieom_base::static_susceptibility_at_high_symmetry(const std::vector<Vector>& states) {
    const auto momenta = square_lattice_high_symmetry_path<L>();
    std::array<std::vector<double>, number_of_high_symmtry_momenta<L>> result;

    for (std::size_t i = 0U; i < number_of_high_symmtry_momenta<L>; ++i) {
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

bool ieom_base::found_pomeranchuk_instability() const noexcept {
    return _found_pomeranchuk;
}

double ieom_base::hartree_fock_filling() const noexcept 
{
    return 2. * std::accumulate(this->occupation_numbers.begin(), this->occupation_numbers.end(), double{}) / N;
}

double ieom_base::chi_number(
    const momentum& k,
    const momentum& kx) const noexcept
{
    const double E  = dispersion[k];
    const double Ex = dispersion[kx];

    const double a = 0.5 * beta * E;
    const double b = 0.5 * beta * Ex;
    const double x = a - b;

    const double log_chi =
        std::log(0.5 * beta)
        + log_sinhc(x)
        - log_cosh(a)
        - log_cosh(b);

    return std::exp(log_chi);
}

double ieom_base::chi_pair(
    const momentum& k,
    const momentum& kx) const noexcept
{
    const double E  = dispersion[k];
    const double Ex = dispersion[kx];

    const double a = 0.5 * beta * E;
    const double b = 0.5 * beta * Ex;
    const double x = a + b;

    const double log_chi =
        std::log(0.25 * beta)
        + log_sinhc(x)
        - log_cosh(a)
        - log_cosh(b);

    return std::exp(log_chi);
}


//-----------------------------------------------------//
// Private methods

double ieom_base::compute_chemical_potential(const Model& model, const double target_filling) {
    constexpr double MU_STEP = 0.1;

    auto fit_occupation = [&model, &target_filling, this](double x) {
        double filling{};
        for (const auto& energy : dispersion) {
            filling += model.fermi_function(energy - x);
        }
        filling *= 2. / N;
        return target_filling - filling;
    };

    double definitly_too_small = chemical_potential;
    double definitly_too_large = chemical_potential;

    double current_deviation = fit_occupation(chemical_potential);
    // In the unlikely case we have a perfect hit, we don't need to work.
    if (is_zero(current_deviation)) return chemical_potential;

    // We are in a non-symmetry-broken state, so the filling is a monotone function of mu 
    if (current_deviation < 0.0) {
        // filling > target_filling => mu must be lowered
        definitly_too_large = chemical_potential;
        do {
            definitly_too_small -= MU_STEP;
            current_deviation = fit_occupation(definitly_too_small);
            if (definitly_too_small < -16.) {
                throw std::runtime_error("definitly_too_small ran below -2 * band width!");
            }
        } while (current_deviation < 0.0);
        if (is_zero(current_deviation)) {
            chemical_potential = definitly_too_small;
            return definitly_too_small;
        }
    }
    else {
        // filling < target_filling => mu must be increased
        definitly_too_small = chemical_potential;
        do {
            definitly_too_large += MU_STEP;
            current_deviation = fit_occupation(definitly_too_large);
            if (definitly_too_large > 16.) {
                throw std::runtime_error("definitly_too_large ran above 2 * band width!");
            }
        } while (current_deviation > 0.0);
        if (is_zero(current_deviation)) {
            chemical_potential = definitly_too_large;
            return definitly_too_large;
        }
    }

    std::uintmax_t boost_max_it{100U};
    // boost::math::tools::eps_tolerance<double>(36) corresponds to about 2.91e-11
    const auto best_mu =
         boost::math::tools::toms748_solve(fit_occupation, definitly_too_small, definitly_too_large,
                                              boost::math::tools::eps_tolerance<double>(36), boost_max_it);

    const double low  = std::abs(fit_occupation(best_mu.first));
    const double mid  = std::abs(fit_occupation(0.5 * (best_mu.first + best_mu.second)));
    const double high = std::abs(fit_occupation(best_mu.second));

    double mu = best_mu.first;
    double min_deviation = low;
            
    if (mid < min_deviation) {
        mu = 0.5 * (best_mu.first + best_mu.second);
        min_deviation = mid;
    }
    
    if (high < min_deviation) {
        mu = best_mu.second;
    }
    
    chemical_potential = mu;
    return mu;
}

void ieom_base::solve_hartree_fock(const Model& model) {
    const auto compute_occupation_numbers = [this, &model](const Vector& dispersion) {
        for (int i = 0; i < N; ++i) {
            occupation_numbers[i] = model.fermi_function(dispersion[i] - chemical_potential);
        }
    };

    Eigen::Map<const Vector> epsilon_tilde_map(epsilon_tilde.data(), N);
    Vector self_energy = Eigen::Map<const Vector>(flow_result.dispersion.data(), N) 
        - epsilon_tilde_map;

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
        for (int i = 0; i < N; ++i) {
            dispersion[i] = epsilon_tilde[i] + self_energy(i);
        }
        chemical_potential = compute_chemical_potential(model, model.filling);
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
    
    for (int i = 0; i < N; ++i) {
        dispersion[i] = epsilon_tilde[i] + self_energy(i);
    }
    compute_occupation_numbers(epsilon_tilde_map + self_energy);
    for (int i = 0; i < N; ++i) {
        dispersion[i] = epsilon_tilde[i] + self_energy(i) - chemical_potential;
    }

    for (momentum k = momentum::begin(); k != momentum::end(); ++k) {
        const momentum kprime = k.rotate_90();
        if (!float_equal(dispersion[k], dispersion[kprime])) {
            _found_pomeranchuk = true;
            break;
        }
    }
}

void ieom_base::prepare_static(const momentum& x) {
    fill_matrices(x);
    const double hermiticity_error = std::sqrt((U_matrix - U_matrix.adjoint()).squaredNorm() / U_matrix.squaredNorm());
    if (is_non_zero(hermiticity_error)) {
        throw std::runtime_error("U_matrix is not Hermitian! ||U-U^+|| / ||U|| = " + std::to_string(hermiticity_error));
    }
    for (auto& chi : chi_matrix) {
        chi = std::sqrt(chi);
    }
    eigen_solver.compute(chi_matrix.asDiagonal() * U_matrix * chi_matrix.asDiagonal());
}

ChiContainer::ChiContainer(const Eigen::VectorXd& eigen_energies, std::size_t n_states)
    : energies(eigen_energies.data(), eigen_energies.data() + eigen_energies.size()) {
    weights.reserve(n_states);
}

void to_json(nlohmann::json& j, const ChiContainer& container) noexcept {
    j = nlohmann::json{{"energies", container.energies}, {"weights", container.weights}};
}

}  // namespace NickelCUT::susceptibility