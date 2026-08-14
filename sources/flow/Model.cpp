#include "Model.hpp"
#include "../helper_functions.hpp"

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
    if (this->beta >= 0.0) {
        return 1. / (1. + std::exp(beta*energy));
    }
    if (is_zero(energy)) {
        return 0.5;
    }
    return energy < 0.0 ? 1.0 : 0.0;
}

} // namespace NickelCUT::flow
