#pragma once

namespace NickelCUT::flow
{

struct Model {
    double U_0;
    double tprime;
    double chemical_potential;
    double beta;

    double epsilon_0(double kx, double ky) const noexcept;

    double fermi_function(double energy) const noexcept;
};

} // namespace NickelCUT::flow
