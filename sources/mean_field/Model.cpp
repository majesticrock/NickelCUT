#include "Model.hpp"

#include "../flow/momentum_iterator.hpp"
#include "../L.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>

#include <iostream>
#include <filesystem>
#include <fstream>

#ifdef NICKEL_COMPLEX
#define __conj(x) std::conj(x)
#else
#define __conj(x) x
#endif

namespace NickelCUT::mean_field
{
    
Model::Model(const std::string& flow_state_dir)
    : flow::Model{0.5, 0.0, 0.0, -1.0}, // TODO: Placeholder values!
    deltas( ModelAttributes<double>::Random(5*N)),
    initial_filling{0.0}
{
    const std::string state_file = flow_state_dir + "lowest_ROD_state.bin";
    if(!std::filesystem::exists(state_file)) {
        throw std::runtime_error("FileNotFound: " + state_file);
    }
    std::cout << "Loading flow state from " << state_file << std::endl;

    {
        std::ifstream ifs(state_file, std::ios::binary);
        if (ifs.good()) {
            boost::archive::binary_iarchive ia(ifs);
            ia >> this->flow_state;
        } else {
            throw std::runtime_error("Inputstream for " + state_file + " is bad!");
        }
    }

    for (auto p = flow::momentum_iterator<L>::begin(); p != flow::momentum_iterator<L>::end(); ++p) {
        initial_filling += 2 * fermi_function(flow_state.dispersion[p]); // factor 2 for spin degeneracy
    }
    initial_filling /= N;
}

void Model::iteration_step(const ParameterVector& initial_values, ParameterVector& result) 
{
    result.setZero();
    this->deltas.fill_with(initial_values);

    // chemical potential, somehow?
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

void Model::fill_hamiltonian(const momentum_t& p)
{
    hamiltonian(0, 0) = flow_state.epsilon_tilde[p] + epsilon_I_up(p) - chemical_potential;
    hamiltonian(0, 1) = Delta_DW_up(p);
    hamiltonian(0, 2) = __conj(Delta_SC(p));
    hamiltonian(0, 3) = 0.0;

    hamiltonian(1, 0) = Delta_DW_up(p + flow::Q<L>);
    hamiltonian(1, 1) = flow_state.epsilon_tilde[p + flow::Q<L>] + epsilon_I_up(p + flow::Q<L>) - chemical_potential;
    hamiltonian(1, 2) = 0.0;
    hamiltonian(1, 3) = __conj(Delta_SC(p + flow::Q<L>));

    hamiltonian(2, 0) = Delta_SC(p);
    hamiltonian(2, 1) = 0.0;
    hamiltonian(2, 2) = -flow_state.epsilon_tilde[-p] + epsilon_I_down(-p) + chemical_potential;
    hamiltonian(2, 3) = -Delta_DW_down(-p);

    hamiltonian(3, 0) = 0.0;
    hamiltonian(3, 1) = Delta_SC(p + flow::Q<L>);
    hamiltonian(3, 2) = -Delta_DW_down(flow::Q<L> - p);
    hamiltonian(3, 3) = -flow_state.epsilon_tilde[flow::Q<L> - p] + epsilon_I_down(flow::Q<L> - p) + chemical_potential;
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

} // namespace NickelCUT::mean_field
