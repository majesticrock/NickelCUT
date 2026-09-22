#include "Model.hpp"

#include "../flow/momentum_iterator.hpp"
#include "../flow/Model.hpp"
#include "../L.hpp"
#include "../flow/flow_state_serialization.hpp"
#include "../flow/data_file_names.hpp"

#include <iostream>
#include <filesystem>
#include <fstream>
#include <set>
#include <limits>

#ifdef NICKEL_COMPLEX
#define __conj(x) std::conj(x)
#else
#define __conj(x) x
#endif

constexpr double MU_TOL = 1e-12;
constexpr double MU_INITIAL_STEP = 0.1;

namespace NickelCUT::mean_field
{
    
Model::Model(const std::string& binary_data_dir, double U_0_, double tprime_, double E_F_, double temperature_)
    : flow::Model(U_0_, tprime_, E_F_, temperature_),
    extracted_channels(flow::deserialize_extracted_channels(binary_data_dir + flow::Model::data_dir_name(), flow::data_file_names::LOWEST_ROD_EXTRACTED_CHANNELS)),
    deltas(5*N, 0.0),
    chemical_potential{0.0}
{
    filling = 0.0;
    // Compute the filling without any mean-field order
    for (std::size_t p=0U; p < extracted_channels.dispersion.size(); ++p) {
        filling += 2 * fermi_function_zero_temperature(extracted_channels.dispersion[p]); // factor 2 for spin degeneracy
    }
    filling /= N;
    std::cout << "Filling = " << filling << std::endl;

    for (momentum_t p = momentum_t::begin(); p != momentum_t::end(); ++p) {
        Sigma_up(p) = extracted_channels.dispersion[p] - extracted_channels.epsilon_tilde[p];
        Sigma_down(p) = Sigma_up(p);
    }
    //compute_chemical_potential();

    for (int p = 0; p < N; ++p) {
        // Some mindlessly picked initial values
        Delta_SC(p) = 1. - (N/2 - p) / 100.;
        Delta_DW_up(p) = 1. - (N/2 - p) / 100.;
        Delta_DW_down(p) = -1. + (N/2 - p) / 100.;
    }
}

void Model::iteration_step(const ParameterVector& initial_values, ParameterVector& result) 
{
    static int i=0; ++i;
    result.setZero();
    this->deltas.fill_with(initial_values);
    
    //compute_chemical_potential();
    // Maybe add the chemical potential to the self-consistency values

    for (momentum_t q = momentum_t::begin(); q != momentum_t::half_end(); ++q) {
        compute_rho(q); // rho(row, col)
        for (int k = 0; k < N; ++k) {
            // Delta_SC
            result(k) -= extracted_channels.superconductivity(k, q) * rho(0, 2)
                + extracted_channels.superconductivity(k, flow::PI<L> + q) * rho(1, 3);
            
            // Delta_DW_up
            result(k + N) += extracted_channels.density_wave.first(k, -q) * rho(2, 3)
                + extracted_channels.density_wave.first(k, flow::PI<L> - q) * rho(3, 2);
            result(k + N) -= extracted_channels.density_wave.second(k, q) * rho(1, 0) 
                + extracted_channels.density_wave.second(k, flow::PI<L> + q) * rho(0, 1);

            // Delta_DW_down
            result(k + 2*N) -= extracted_channels.density_wave.first(k, q) * rho(1, 0)
                + extracted_channels.density_wave.first(k, flow::PI<L> + q) * rho(0, 1);
            result(k + 2*N) += extracted_channels.density_wave.second(k, -q) * rho(2, 3) 
                + extracted_channels.density_wave.second(k, flow::PI<L> - q) * rho(3, 2);

            // Sigma_up
            result(k + 3*N) += extracted_channels.single_particle_energy.first(k, -q) * rho(2, 2)
                + extracted_channels.single_particle_energy.first(k, flow::PI<L> - q) * rho(3, 3);
            result(k + 3*N) += extracted_channels.single_particle_energy.second(k, q) * (1. - rho(0, 0)) 
                + extracted_channels.single_particle_energy.second(k, flow::PI<L> + q) * (1. - rho(1, 1));

            // Sigma_down
            result(k + 4*N) += extracted_channels.single_particle_energy.first(k, q) * (1. - rho(0, 0))
                + extracted_channels.single_particle_energy.first(k, flow::PI<L> + q) * (1. - rho(1, 1));
            result(k + 4*N) += extracted_channels.single_particle_energy.second(k, -q) * rho(2, 2) 
                + extracted_channels.single_particle_energy.second(k, flow::PI<L> - q) * rho(3, 3);
        }
    }

    this->deltas.fill_with(result, 0.5);
    this->deltas.clear_noise(PRECISION);

    std::cout << i << ": " << max_Delta_SC() << "\t" << max_Delta_AFM() << "\t" << max_Delta_CDW() 
        << "\t" << chemical_potential << "\t" << compute_filling() << std::endl;

    result -= initial_values;
}

double Model::compute_filling()
{
    double filling{};
    for (momentum_t p = momentum_t::begin(); p != momentum_t::half_end(); ++p) {
        compute_rho(p);
        filling += 2 - rho(0,0) - rho(1,1) + rho(2,2) + rho(3,3);
    }
    filling /= N;
    return filling;
}

double Model::max_Delta_SC() const noexcept 
{
    double current{};
    for (std::size_t i=0U; i<N; ++i) {
        if (std::abs(Delta_SC(i)) > current) 
            current = std::abs(Delta_SC(i));
    }
    return current;
}

double Model::max_Delta_AFM() const noexcept {
    double current{};
    for (std::size_t i=0U; i<N; ++i) {
        if (std::abs(Delta_DW_up(i) - Delta_DW_down(i)) > current) 
            current = std::abs(Delta_DW_up(i) - Delta_DW_down(i));
    }
    return 0.5 * current;
}

double Model::max_Delta_CDW() const noexcept {
    double current{};
    for (std::size_t i=0U; i<N; ++i) {
        if (std::abs(Delta_DW_up(i) + Delta_DW_down(i)) > current) 
            current = std::abs(Delta_DW_up(i) + Delta_DW_down(i));
    }
    return 0.5 * current;
}

void Model::fill_hamiltonian(const momentum_t& p)
{
    // h(row, col)
    hamiltonian(0, 0) = dispersion_up(p);
    hamiltonian(0, 1) = Delta_DW_up(p);
    hamiltonian(0, 2) = __conj(Delta_SC(p));
    hamiltonian(0, 3) = 0.0;

    hamiltonian(1, 0) = Delta_DW_up(p + flow::PI<L>);
    hamiltonian(1, 1) = dispersion_up(p + flow::PI<L>);
    hamiltonian(1, 2) = 0.0;
    hamiltonian(1, 3) = __conj(Delta_SC(p + flow::PI<L>));

    hamiltonian(2, 0) = Delta_SC(p);
    hamiltonian(2, 1) = 0.0;
    hamiltonian(2, 2) = -dispersion_down(-p);
    hamiltonian(2, 3) = -Delta_DW_down(-p);

    hamiltonian(3, 0) = 0.0;
    hamiltonian(3, 1) = Delta_SC(p + flow::PI<L>);
    hamiltonian(3, 2) = -Delta_DW_down(flow::PI<L> - p);
    hamiltonian(3, 3) = -dispersion_down(flow::PI<L> - p);
}

void Model::compute_rho(const momentum_t& p) 
{
    fill_hamiltonian(p);
    eigensolver.compute(hamiltonian);
    rho = eigensolver.eigenvectors()
        * eigensolver.eigenvalues().unaryExpr([this](const double E) {
            return 1. - this->fermi_function(E);
        }).asDiagonal()
        * eigensolver.eigenvectors().adjoint();
}

void Model::compute_chemical_potential()
{
    //static int i=0; ++i;
    // If the chemical potential is too small, the computed filling will be too small
    // and therefore this lambda will return -1
    // If the chemical potential is too large, the lambda will return 1
    // If the chemical potential is just right, the lambda returns 0
    auto filling_func = [this](double mu_) {
        this->chemical_potential = mu_;
        return compute_filling() - filling;
    };
    double low{ chemical_potential - MU_INITIAL_STEP }; 
    double up{ chemical_potential + MU_INITIAL_STEP };

    std::set<double> unique_energies;
    for(momentum_t p = momentum_t::begin(); p != momentum_t::half_end(); ++p) {
        fill_hamiltonian(p);
        eigensolver.compute(hamiltonian, Eigen::EigenvaluesOnly);
        for (const auto& ev : eigensolver.eigenvalues()) {
            unique_energies.insert(ev);
        }
    }
    for (const auto& energy : unique_energies) {
        if (is_zero(filling_func(energy))) {
            // The chemical potential should lie exactly on one of the energy levels
            return;
        }
    }
    // The chemical potential should lie somewhere between the energy levels

    // Starting values are for now guess work
    double f_low{filling_func(low)};
    while (greater_or_almost_equal(f_low, 0.0)) {
        low -= MU_INITIAL_STEP;
        f_low = filling_func(low);
    }

    double f_up{filling_func(up)};
    while (!greater_or_almost_equal(f_up, 0.0)) {
        up += MU_INITIAL_STEP;
        f_up = filling_func(up);
    }
    const double initial_up{up};

    double center{ 0.5 * (up + low) };
    double f_center{filling_func(center)};

    // First try to find mu_-; the smallest value for the chemical potential that matches our desired filling
    while(up - low > MU_TOL) {
        if (greater_or_almost_equal(f_center, 0.0)){
            up = center;
            f_up = f_center;
        }
        else {
            low = center;
            f_low = f_center;
        }
        center = 0.5 * (up + low);
        f_center = filling_func(center);
    }
    const double mu_minus = center;

    // Now try to find mu_+; the largest value for the chemical potential that matches our desired filling
    low = mu_minus;
    up = initial_up;
    while (less_or_almost_equal(f_up, 0.0)) {
        up += MU_INITIAL_STEP;
        f_up = filling_func(up);
    }
    center = 0.5 * (up + low);

    f_low = filling_func(low);
    f_up = filling_func(up);
    f_center = filling_func(center);
    while(up - low > MU_TOL) {
        if (less_or_almost_equal(f_center, 0.0)) {
            low = center;
            f_low = f_center;
        }
        else {
            up = center;
            f_up = f_center;
        }
        center = 0.5 * (up + low);
        f_center = filling_func(center);
    }
    const double mu_plus = center;

    this->chemical_potential = 0.5 * (mu_minus + mu_plus);
    const double reached_filling = compute_filling();
    if (!is_zero(reached_filling - filling)) {
        throw std::runtime_error("Failed finding the chemical potential! Wanted a filling of " + std::to_string(filling) + " but got " + std::to_string(reached_filling));
    }
}

} // namespace NickelCUT::mean_field
