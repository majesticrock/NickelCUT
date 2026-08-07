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
    coeff_array occupation_numbers;
    double constant{};

    FlowContainer() : constant{} 
    {
        for (int x=0; x<L; ++x) {
            const double kx = std::numbers::pi * (static_cast<double>(x) / static_cast<double>(L) - 1.0);
            for (int y=0; y<L; ++y) {
                const double ky = std::numbers::pi * (static_cast<double>(y) / static_cast<double>(L) - 1.0);
                dispersion[x*L + y] = epsilon_0(kx, ky);
                occupation_numbers[x*L + y] = fermi_function(dispersion[x*L + y]);
            }
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