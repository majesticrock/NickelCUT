#include "Verifier.hpp"

#include "../commutator_of_cut.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>

using namespace mrock::symbolic_operators;
using namespace mrock::symbolic_operators::experimental;

namespace NickelCUT::commute::verify
{

Verifier::Verifier()
    : occupations(make_default_occupations()),
     cosines(make_cosines()),
     templates(NickelCUT::commute::get_wick_templates()),
     symmetries(NickelCUT::commute::get_symmetries()),
     identity(make_identity())
{
    momenta.fill(IntMomentum<L>(0));
    indices.fill(0);
    indices[static_cast<unsigned char>(Index::SpinUp)] = 0;
    indices[static_cast<unsigned char>(Index::SpinDown)] = 1;
}

Verifier::Verifier(const std::array<double, L>& _occupations)
    : occupations(_occupations),
     cosines(make_cosines()),
     templates(NickelCUT::commute::get_wick_templates()),
     symmetries(NickelCUT::commute::get_symmetries()),
     identity(make_identity())
{
    momenta.fill(IntMomentum<L>(0));
    indices.fill(0);
    indices[static_cast<unsigned char>(Index::SpinUp)] = 0;
    indices[static_cast<unsigned char>(Index::SpinDown)] = 1;
}

Verifier::SparseMatrix Verifier::symbolic_to_matrix(const TermCollector& terms) const
{
    return symbolic_to_matrix_impl(terms);
}

Verifier::SparseMatrix Verifier::symbolic_to_matrix(const WickOrderedCollector& terms) const
{
    return symbolic_to_matrix_impl(terms);
}

bool Verifier::operator()(const TermCollector& expression, const char* name, bool hermitian) const 
{
    SparseMatrix before;
    return operator()(expression, name, hermitian, before);
}

bool Verifier::operator()(const TermCollector& expression, const char* name, bool hermitian, SparseMatrix& save_matrix) const {
    save_matrix = symbolic_to_matrix(expression);
    if (hermitian) {
        const SparseMatrix transp = save_matrix.transpose();
        const double diff = (save_matrix - transp).norm();
        if (diff > 1e-12) {
            std::cerr << name << " is not Hermitian, but it should be! Diff = " << diff << "\n";
            std::cerr << "\\begin{align*}\n" << expression << "\\end{align*}" << std::endl;
            return false;
        }
    }
    WickOrderedCollector after = wick_decompose(expression, templates);
    clean_wick_ordered_terms(after, symmetries);
    
    const SparseMatrix after_matrix = symbolic_to_matrix(after);
    if (!matrices_equal(save_matrix, after_matrix)) {
        std::cerr << name << " changed under Wick decomposition. Diff = " << (save_matrix - after_matrix).norm() << "\n";
        std::cerr << "\\begin{align*}\n" << expression << "\\end{align*}\n\\begin{align*}" << after << "\\end{align*}" << std::endl;  
        return false;
    }
    return true;
}

bool Verifier::operator()(const WickOrderedCollector& expression, const char* name, const SparseMatrix compare) const {
    const SparseMatrix expression_matrix = symbolic_to_matrix(expression);
    if (!matrices_equal(compare, expression_matrix)) {
        std::cerr << name << " changed. Diff = " << (compare - expression_matrix).norm() << "\n";
        std::cerr << "\\begin{align*}\n" << expression << "\\end{align*}" << std::endl;  
        return false;
    }
    return true;
}

std::array<double, L> Verifier::make_cosines() {
    std::array<double, L> ret{};
    for (int i=0; i<L; ++i) {
        ret[i] = std::cos(std::numbers::pi * ((2. * i / L) - 1.));
    }
    return ret;
}

std::array<double, L> Verifier::make_default_occupations() {
    std::array<double, L> ret{};
    for (int i=0; i<L; ++i) {
        ret[i] = std::cos(std::numbers::pi * (2. * i / L - 1.));
        ret[i] = std::abs(ret[i]) < 1e-12 ? 0.5 :
                    (ret[i] < 0. ? 1. : 0.);
    }
    return ret;
}

Verifier::SparseMatrix Verifier::make_identity() {
    SparseMatrix ret(matrix_size, matrix_size);
    ret.setIdentity();
    return ret;
}

bool Verifier::matrices_equal(const SparseMatrix& lhs, const SparseMatrix& rhs) const {
    return (lhs - rhs).norm() <= 1.e-12;
}

Verifier::SparseMatrix Verifier::operator_string(const std::vector<int>& modes, const std::vector<bool>& creation) const 
{
    if (modes.size() != creation.size()) {
        throw std::invalid_argument("operator and mode counts differ");
    }
    if (modes.empty()) return identity;

    SparseMatrix result(matrix_size, matrix_size);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(matrix_size);

    for (uint32_t input = 0; input < static_cast<uint32_t>(matrix_size); ++input) {
        uint32_t output = input;
        double sign = 1.;
        bool valid = true;
        for (std::size_t i = modes.size(); i-- > 0;) {
            const uint32_t mask = 1u << modes[i];
            const bool occupied = (output & mask) != 0u;
            if (occupied == creation[i]) {
                valid = false;
                break;
            }
            if (creation[i]) {
                output |= mask;
            } else {
                output &= ~mask;
            }
            if ((std::popcount(output & (mask - 1u)) & 1) != 0) {
                sign = -sign;
            }
        }
        if (valid) {
            triplets.emplace_back(static_cast<int>(output), static_cast<int>(input), sign);
        }
    }
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

Verifier::SparseMatrix Verifier::normal_ordered_operator_string(const std::vector<int>& modes, const std::vector<bool>& creation) const 
{
    if (modes.size() != creation.size()) {
        throw std::invalid_argument("operator and mode counts differ");
    }

    const auto contraction = [&](int l, int r)
    {
        if (modes[l] != modes[r]) return 0.0;
        if (!creation[l] || creation[r]) return 0.0;

        return occupations[modes[l] / 2];
    };

    const auto reduced_string = [&](const std::vector<int>& indices) {
        std::vector<int> r_modes(indices.size());
        std::vector<bool> r_creation(indices.size());
        for (std::size_t i=0U; i<indices.size(); ++i) {
            r_modes[i]    = modes[indices[i]];
            r_creation[i] = creation[indices[i]];
        }
        return normal_ordered_operator_string(r_modes, r_creation);
    };

    const auto all_indices_except = [&](std::set<int>&& skip) {
        std::vector<int> indices(modes.size() - skip.size());
        int c=0;
        for (auto& index : indices) {
            while (std::find(skip.begin(), skip.end(), c) != skip.end()) {
                ++c;
            }
            index = c++;
        }
        return indices;
    };

    if (modes.empty()) {
        return identity;
    } 
    if (modes.size() == 2U) {
        return (operator_string(modes, creation) - contraction(0, 1) * identity);
    }
    if (modes.size() == 4U) {
        SparseMatrix ret = operator_string(modes, creation);

        // Single contractions
        for (int i=0; i<3; ++i) {
            for (int j=i+1; j<4; ++j) {
                // Gives +1 if i and j are adjacent or have an even number of operators between them
                // For instance 0,1 or 0,3
                const double sign = ((j-i) & 1) ? 1.0 : -1.0;
                ret -= sign * contraction(i, j) * reduced_string(all_indices_except({i, j}));
            }
        }

        // double contractions
        double constant=0.;
        for (int j=1; j<4; ++j) {
            const double sign = ((j-1) & 1) ? 1.0 : -1.0;
            const auto remaining_indices = all_indices_except({0, j});
            constant += sign * contraction(0, j) * contraction(remaining_indices[0], remaining_indices[1]);
        }
        ret += constant * identity;

        return ret;
    }
    if (modes.size() == 6U) {
        SparseMatrix ret = operator_string(modes, creation);

        // single contractions
        for (int i = 0; i < 5; ++i) {
            for (int j = i + 1; j < 6; ++j) {
                // Gives +1 if i and j are adjacent or have an even number of operators between them
                // For instance 0,1 or 0,3
                const double sign = ((j-i) & 1) ? 1.0 : -1.0;
                ret -= sign * contraction(i, j) * reduced_string(all_indices_except({i, j}));
            }
        }

        // double contractions
        for (int i = 0; i < 5; ++i) {
            for (int j = i + 1; j < 6; ++j) {
                // i, j select the normal ordered operator string.
                // The remaing 4 indizes are going to be contracted.
                const auto remaining = all_indices_except({i, j});
                const double sign1 = ((j-i) & 1) ? 1.0 : -1.0;

                // Contract the first remaining operator with one of the other three.
                for (int k = 1; k < 4; ++k) {
                    const double sign2 = (k & 1) ? 1.0 : -1.0;
                    const auto inner_remaining = all_indices_except({i, j, remaining[0], remaining[k]});

                    ret -= sign1 * sign2
                         * contraction(inner_remaining[0], inner_remaining[1])
                         * contraction(remaining[0], remaining[k])
                         * reduced_string({i, j});
                }
            }
        }

        // triple contractions
        double constant = 0.0;
        for (int j = 1; j < 6; ++j) {
            const auto remaining = all_indices_except({0, j});
            const double sign1 = (j & 1) ? 1.0 : -1.0;

            for (int k = 1; k < 4; ++k) {
                const double sign2 = (k & 1) ? 1.0 : -1.0;
                const auto inner_remaining = all_indices_except({0, j, remaining[0], remaining[k]});

                constant -= sign1 * sign2 
                          * contraction(0, j)
                          * contraction(remaining[0], remaining[k])
                          * contraction(inner_remaining[0], inner_remaining[1]);
            }
        }
        ret += constant * identity;
        return ret;
    }
    throw std::invalid_argument("Number of operators not implemented. Expected only 0, 2, 4, or 6, but got " + std::to_string(modes.size()) + ".");
}

IntMomentum<L> Verifier::momentum_lookup(const Momentum& momentum) const noexcept {
    IntMomentum<L> value(momentum.add_PI ? 0 : L/2);
    for (const auto& symbol : momentum.momentum_list) {
        value += symbol.factor * momenta[static_cast<unsigned char>(symbol.name)];
    }
    return value;
}

int Verifier::index_lookup(const Index& index) const noexcept {
    if (index == Index::AntiSigma) {
        // AntiSigma is meant to be the anti-parallel spin state to Sigma
        // So if Sigma is spin up (represented by 0), AntiSigma is spin down (represented by 1), and vice versa
        return (int)(!indices[static_cast<unsigned char>(Index::Sigma)]);
    }
    if (index == Index::AntiSigmaPrime) {
        return (int)(!indices[static_cast<unsigned char>(Index::SigmaPrime)]);
    }
    return indices[static_cast<unsigned char>(index)];
}

double Verifier::coefficient_value(const Coefficient& coefficient) const
{
    const IntMomentum<L> k = momentum_lookup(coefficient.momenta[0]);
    if (coefficient.name == "\\tilde{\\varepsilon}") {
        return cosines[k];
    }
    const IntMomentum<L> p = momentum_lookup(coefficient.momenta[1]);
    const IntMomentum<L> q = momentum_lookup(coefficient.momenta[2]);
    
    double interaction = cosines[k] * cosines[p] + cosines[k+q] * cosines[p-q];
    if (coefficient.indices.size() == 2U) {
        interaction *= (index_lookup(coefficient.indices[0]) == index_lookup(coefficient.indices[1]) ? 1. : cosines[q]);
    }
    else if (coefficient.indices.size() == 1U) {
        interaction *= (coefficient.indices[0] == Index::Parallel ? 1. : cosines[q]);
    }
    if (coefficient.name == "U") {
        return interaction;
    }
    if (coefficient.name == "\\alpha") {
        return interaction * (cosines[k] + cosines[p] - cosines[k + q] - cosines[p - q]);
    }
    throw std::runtime_error("unknown coefficient: " + coefficient.name);
}

} // NickelCUT::commute::verify