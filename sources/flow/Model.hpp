#pragma once

#include <nlohmann/json_fwd.hpp>

#include <string>

namespace NickelCUT::flow
{

struct Model {
    const double U_0{};
    const double tprime{};
    const double E_F{};
    const double temperature{};
    
    const double beta{};
    //double chemical_potential{};
    double filling{};

    Model(double U_0_, double tprime_, double E_F_, double temperature_);

    double epsilon_0(double kx, double ky) const noexcept;

    double fermi_function_zero_temperature(double energy) const noexcept;

    double fermi_function(double energy) const noexcept;

    std::string info() const;

    std::string data_dir_name() const;

    nlohmann::json generate_meta_data_json() const;
};

} // namespace NickelCUT::flow
