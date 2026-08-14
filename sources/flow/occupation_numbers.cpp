#include "occupation_numbers.hpp"

#include "momentum_iterator.hpp"
#include "Model.hpp"
#include "FlowContainer.hpp"

namespace NickelCUT::flow {

std::array<double, N> occupation_numbers;

void compute_occupation_numbers(const Model& model, const FlowContainer& init) {
    for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
        occupation_numbers[p.get_position()] = model.fermi_function(init.dispersion[p.get_position()]);
    }
}

}