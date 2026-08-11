#include "Model.hpp"

#include <cmath>

namespace NickelCUT::flow
{
double Model::epsilon_0(double kx, double ky) const noexcept
{
    return -2. * (
                std::cos(kx) + std::cos(ky) 
                + 2. * tprime * std::cos(kx) * std::cos(ky) 
            ) - chemical_potential;
}

double Model::fermi_function(double energy) const noexcept
{
    return 1. / (1. + std::exp(beta*energy));
}

} // namespace NickelCUT::flow
