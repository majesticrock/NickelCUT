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

double Model::fermi_function_zero_temperature(double energy) const noexcept 
{
    if (is_zero(energy)) {
        return 0.5;
    }
    return energy < 0.0 ? 1.0 : 0.0;
}

double Model::fermi_function(double energy) const noexcept
{
    if (this->beta >= 0.0) {
        return 1. / (1. + std::exp(beta*energy));
    }
    return fermi_function_zero_temperature(energy);
}

std::string Model::info() const noexcept
{
    return "NickelCUT Model: U_0=" + std::to_string(U_0)
        + "   tprime=" + std::to_string(tprime)
        + "   chemical_potential=" + std::to_string(chemical_potential)
        + "   beta=" + std::to_string(beta); 
}

} // namespace NickelCUT::flow
