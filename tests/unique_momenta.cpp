#include "../sources/flow/unique_momenta.hpp"

#include <array>
#include <cmath>
#include <iostream>

using namespace NickelCUT::flow;

int main() {
    constexpr int L=20;
    std::array<double, L*L> epsilon;
    for (momentum_iterator<L> k=momentum_iterator<L>::begin(); k != momentum_iterator<L>::end(); ++k) {
        epsilon[k] = -2. * (std::cos(k.get_kx()) + std::cos(k.get_ky()))
            + 0.5 * std::cos(k.get_kx()) * std::cos(k.get_ky())
            - 3.;
    }

    for (momentum_iterator<L> k=momentum_iterator<L>::begin(); k != momentum_iterator<L>::end(); ++k) {
        if (std::abs(epsilon[k] - epsilon[unique_momenta<L>::reduce(k)]) > 1e-15) {
            std::cerr << "epsilon_k changed under the transformation!\n"
                << "eps_k = " << epsilon[k] << "     eps_k' = " << epsilon[unique_momenta<L>::reduce(k)] << "\n"
                << "k = " << k << "    k' = " << momentum_iterator<L>(unique_momenta<L>::reduce(k))
                << std::endl;
            return 1;
        }
    }

    for (const auto& u_it : unique_momenta<L>::uniques) {
        if (u_it != momentum_iterator<L>(unique_momenta<L>::reduce(u_it))) {
            std::cerr << "unqiue k changed under the transformation!\n"
                << "k = " << u_it << "    k' = " << momentum_iterator<L>(unique_momenta<L>::reduce(u_it))
                << std::endl;
        }
    }

    for (momentum_iterator<L> k=momentum_iterator<L>::begin(); k != momentum_iterator<L>::end(); ++k) {
        const int pos_transformed = unique_momenta<L>::reduce(k);
        const momentum_iterator<L> transformed(pos_transformed);
        if (std::find(unique_momenta<L>::uniques.begin(), unique_momenta<L>::uniques.end(), transformed) 
            == unique_momenta<L>::uniques.end()) {
            std::cerr << "k' is not in the unique momenta array!\n"
                << "k = " << k << "    k' = " << transformed << "   " << pos_transformed
                << std::endl;
            return 1;
        }
    }

    momentum_iterator<L> test(0);
    if (test._pos != unique_momenta<L>::reduce(test)) {
        std::cerr << "PI changed under the transformation!\n"
            << "k = " << test << "     k' = " << momentum_iterator<L>(unique_momenta<L>::reduce(test)) << std::endl;
        return 1;
    }

    test = momentum_iterator<L>(L/2, L/2);
    if (test._pos != unique_momenta<L>::reduce(test)) {
        std::cerr << "Gamma changed under the transformation!\n"
            << "k = " << test << "     k' = " << momentum_iterator<L>(unique_momenta<L>::reduce(test)) << std::endl;
        return 1;
    }

    test = momentum_iterator<L>(1, L/2);
    if (momentum_iterator<L>(L/2, 1)._pos != unique_momenta<L>::reduce(test)) {
        std::cerr << "(1,L/2) changed under the transformation!\n"
            << "k = " << test << "     k' = " << momentum_iterator<L>(unique_momenta<L>::reduce(test)) << std::endl;
        return 1;
    }

    return 0;
}