#include "../sources/flow/momentum_iterator.hpp"

#include <array>
#include <cmath>
#include <numbers>

constexpr double U_0 = 0.5;
constexpr double tprime = 0.5;
constexpr double chemical_potential = 0.;
constexpr double beta = 10;

constexpr int L = 6;
constexpr int N = L*L;

using coeff_array = std::array<double, N>;
using mom_it = NickelCUT::flow::momentum_iterator<L>;

double epsilon_0(double kx, double ky) {
    return -2. * (
                std::cos(kx) + std::cos(ky) 
                + 4. * tprime * std::cos(kx) * std::cos(ky) 
            ) - chemical_potential;
}

double fermi_function(double eps) {
    return 1. / (std::exp(beta*eps) + 1);
}

struct FlowContainer {
    std::array<std::array<coeff_array, N>, N> interactions_same_spin;
    std::array<std::array<coeff_array, N>, N> interactions_differing_spin;
    coeff_array dispersion;
    coeff_array epsilon_tilde; // TODO!
    coeff_array occupation_numbers;
    double constant{};

    FlowContainer() : constant{} 
    {
        for (mom_it it = mom_it::begin(); it != mom_it::end(); ++it) {
            dispersion[it.get_position()] = epsilon_0(it.get_kx(), it.get_ky());
            occupation_numbers[it.get_position()] = fermi_function(dispersion[it.get_position()]);
        }

        for (int i=0; i<N; ++i) {
            for (int j=0; j<N; ++j) {
                for (int k=0; k<N; ++k) {
                    interactions_same_spin[i][j][k] = U_0;
                    interactions_differing_spin[i][j][k] = 0.;
                }
            }
        }
    }
};

void flow_step(FlowContainer& result, const FlowContainer& current) {
// Generate this via commute?
}

int main() {

}