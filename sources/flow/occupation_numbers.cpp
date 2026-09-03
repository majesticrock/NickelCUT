#include "occupation_numbers.hpp"

#include "momentum_iterator.hpp"
#include "Model.hpp"

namespace NickelCUT::flow {

std::array<double, N> occupation_numbers;

double compute_occupation_numbers(const Model& model) {
    // model.epsilon_0(p) includes the net chemical potential originating from a finite U
    // Note that the Fock contribution at l=0 vanishes because U_parallel (l=0) = 0.
    // The filling computed via these occupation numbers is therefore the same as the filling
    // of the system with U=0.
    
    double filling{};
    for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
        occupation_numbers[p.get_position()] = model.fermi_function(model.epsilon_0(p.get_kx(), p.get_ky()));
        filling += occupation_numbers[p.get_position()];
    }
    filling /= N;
    return filling;
}

}