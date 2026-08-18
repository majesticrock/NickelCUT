#pragma once

#include <string>

namespace NickelCUT::flow
{

struct Model {
    const double U_0{};
    const double tprime{};
    double chemical_potential{};
    const double beta{};

    double epsilon_0(double kx, double ky) const noexcept;

    double fermi_function_zero_temperature(double energy) const noexcept;

    double fermi_function(double energy) const noexcept;

    std::string info() const noexcept;
};

} // namespace NickelCUT::flow
