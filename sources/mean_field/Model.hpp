#pragma once

#include "../flow/Model.hpp"
#include "../flow/ExtractionContainer.hpp"
#include "../flow/momentum_iterator.hpp"
#include "ModelAttributes.hpp"
#include "../L.hpp"

#include <Eigen/Dense>

#include <cassert>
#include <string>

namespace NickelCUT::mean_field {

struct Model : public flow::Model {
    typedef Eigen::VectorXd ParameterVector;
    typedef flow::momentum_iterator<L> momentum_t;

    flow::ExtractionContainer extracted_channels;
    ModelAttributes<double> deltas;

    double chemical_potential;

    Model(const std::string& binary_data_dir, double U_0_, double tprime_, double E_F_, double temperature_);

    // Functionality for the self-consistency procedure
    void iteration_step(const ParameterVector& initial_values, ParameterVector& result);
    double compute_filling();

    double max_Delta_SC() const noexcept;
    double max_Delta_AFM() const noexcept;
    double max_Delta_CDW() const noexcept;

    // Accessor convenience functions
    inline double Delta_SC(std::size_t i) const noexcept {
        assert(i < N); return deltas[i];
    }
    inline double Delta_DW_up(std::size_t i) const noexcept {
        assert(i < N); return deltas[i + N];
    }
    inline double Delta_DW_down(std::size_t i) const noexcept {
        assert(i < N); return deltas[i + 2*N];
    }
    inline double Sigma_up(std::size_t i) const noexcept {
        assert(i < N); return deltas[i + 3*N];
    }
    inline double Sigma_down(std::size_t i) const noexcept {
        assert(i < N); return deltas[i + 4*N];
    }
    // return reference
    inline double& Delta_SC(std::size_t i) noexcept {
        assert(i < N); return deltas[i];
    }
    inline double& Delta_DW_up(std::size_t i) noexcept {
        assert(i < N); return deltas[i + N];
    }
    inline double& Delta_DW_down(std::size_t i) noexcept {
        assert(i < N); return deltas[i + 2*N];
    }
    inline double& Sigma_up(std::size_t i) noexcept {
        assert(i < N); return deltas[i + 3*N];
    }
    inline double& Sigma_down(std::size_t i) noexcept {
        assert(i < N); return deltas[i + 4*N];
    }
    
    inline double dispersion_up(std::size_t p) const noexcept {
        return extracted_channels.epsilon_tilde[p] + Sigma_up(p) - chemical_potential;
    }
    inline double dispersion_down(std::size_t p) const noexcept {
        return extracted_channels.epsilon_tilde[p] + Sigma_down(p) - chemical_potential;
    }
private:  
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigensolver;
    Eigen::Matrix4d hamiltonian;
    Eigen::Matrix4d rho;

    void fill_hamiltonian(const momentum_t& p);
    void compute_rho(const momentum_t& p);
    void compute_chemical_potential();
};

} // namespace NickelCUT::mean_field
