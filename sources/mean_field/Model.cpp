#include "Model.hpp"

#include "../flow/momentum_iterator.hpp"
#include "../flow/Model.hpp"
#include "../L.hpp"
#include "../flow/flow_state_serialization.hpp"
#include "../flow/data_file_names.hpp"
#include "../helper_functions.hpp"

#include <mrock/utility/InputFileReader.hpp>

#include <algorithm>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <vector>

#ifdef NICKEL_COMPLEX
#define __conj(x) std::conj(x)
#else
#define __conj(x) x
#endif

constexpr double MU_TOL = 1e-12;
constexpr double MU_INITIAL_STEP = 0.1;

namespace NickelCUT::mean_field
{
    
Model::Model(const std::string& binary_data_dir, mrock::utility::InputFileReader& input)
    : flow::Model(input),
    extracted_channels(flow::deserialize_extracted_channels(binary_data_dir + flow::Model::data_dir_name(), flow::data_file_names::LOWEST_ROD_EXTRACTED_CHANNELS)),
    deltas(5*N, 0.0),
    chemical_potential{0.0},
    target_filling{2. * filling} // filling is computed by the parent model.
{
    reset_self_consistency_values();
    std::cout << "Target filling = " << target_filling << "     Filling = " << filling << std::endl;
}

void Model::reset_self_consistency_values() noexcept
{
    for (momentum_t p = momentum_t::begin(); p != momentum_t::end(); ++p) {
        Sigma_up(p) = extracted_channels.dispersion[p] - extracted_channels.epsilon_tilde[p];
        Sigma_down(p) = Sigma_up(p);
    }
    compute_chemical_potential();
    filling = compute_filling();

    for (int p = 0; p < N; ++p) {
        // Some mindlessly picked initial values
        Delta_SC(p) = 1. - (N/2 - p) / 100.;
        Delta_DW_up(p) = 1. - (N/2 - p) / 100.;
        Delta_DW_down(p) = -1. + (N/2 - p) / 100.;
    }
}

void Model::iteration_step(const ParameterVector& initial_values, ParameterVector& result) 
{
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

    result -= initial_values;

    //static int i=0; ++i;
    //std::cout << i << ": " << result.norm() << "  ||  " << max_Delta_SC() << "\t" << max_Delta_AFM() << "\t" << max_Delta_CDW() 
    //    << "\t" << chemical_potential << "\t" << compute_filling() << std::endl;
}

double Model::compute_filling()
{
    double _filling{};
    for (momentum_t p = momentum_t::begin(); p != momentum_t::half_end(); ++p) {
        compute_rho(p);
        _filling += 2 - rho(0,0) - rho(1,1) + rho(2,2) + rho(3,3);
    }
    _filling /= N;
    return _filling;
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

OrderType Model::order_type(double mean_field_precision) const noexcept {
    OrderType ret = OrderType::NormalMetal;
    if (max_Delta_SC() > 10. * mean_field_precision) {
        // SC order exists
        double c_s{};
        double c_d{};
        for (momentum_t k = momentum_t::begin(); k != momentum_t::end(); ++k) {
            c_s += Delta_SC(k);
            c_d += Delta_SC(k) * (std::cos(k.get_kx()) - std::cos(k.get_ky()));
        }
        c_s /= N;
        c_d /= N;
        
        if (std::abs(c_s) > 10. * mean_field_precision) ret = ret | OrderType::swave_SC;
        if (std::abs(c_d) > 10. * mean_field_precision) ret = ret | OrderType::dwave_SC;
    }
    if (max_Delta_AFM() > 10. * mean_field_precision) {
        // AFM order exists
        ret = ret | OrderType::AFM;
    }
    if (max_Delta_CDW() > 10. * mean_field_precision) {
        // CDW order exists
        ret = ret | OrderType::CDW;
    }

    return ret;
}

nlohmann::json Model::selfconsistency_to_json() const noexcept {
    std::array<double, N> _Delta_SC; 
    std::array<double, N> _Delta_AFM;
    std::array<double, N> _Delta_CDW;
    std::array<double, N> _epsilon_up;
    std::array<double, N> _epsilon_down;

    for (std::size_t i = 0U; i < _Delta_SC.size(); ++i) {
        _Delta_SC[i] = Delta_SC(i);
        _Delta_AFM[i] = 0.5 * (Delta_DW_up(i) - Delta_DW_down(i));
        _Delta_CDW[i] = 0.5 * (Delta_DW_up(i) + Delta_DW_down(i));
        _epsilon_up[i] = dispersion_up(i);
        _epsilon_down[i] = dispersion_down(i);
    }

    return nlohmann::json{
        { "Delta_SC",     _Delta_SC     },
        { "Delta_AFM",    _Delta_AFM    },
        { "Delta_CDW",    _Delta_CDW    },
        { "epsilon_up",   _epsilon_up   },
        { "epsilon_down", _epsilon_down }
    };
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
    std::vector<double> energy_levels;
    for (momentum_t p = momentum_t::begin(); p != momentum_t::end(); ++p) {
        energy_levels.push_back(dispersion_up(p) + chemical_potential);
        energy_levels.push_back(dispersion_down(p) + chemical_potential);
    }

    if (energy_levels.empty()) {
        return;
    }

    auto sort_and_unique = [](std::vector<double>& values) {
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end(), float_equal), values.end());
    };

    sort_and_unique(energy_levels);

    std::vector<double> chemical_potentials = energy_levels;
    for (auto lower = energy_levels.begin(), upper = std::next(lower);
         upper != energy_levels.end(); ++lower, ++upper) {
        chemical_potentials.push_back(0.5 * (*lower + *upper));
    }
    sort_and_unique(chemical_potentials);

    double best_chemical_potential = chemical_potentials.front();
    double best_filling_difference = std::numeric_limits<double>::max();

    for (const double candidate : chemical_potentials) {
        chemical_potential = candidate;
        const double filling_difference = std::abs(compute_filling() - target_filling);
        if (filling_difference < best_filling_difference) {
            best_chemical_potential = candidate;
            best_filling_difference = filling_difference;
        }
    }

    chemical_potential = best_chemical_potential;
    filling = compute_filling();
}

} // namespace NickelCUT::mean_field
