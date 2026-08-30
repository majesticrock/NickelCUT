#include "../sources/experimental/WickOrderedCollector.hpp"
#include "../sources/commute/commutator_of_cut.hpp"

#include <mrock/symbolic_operators/Commutation>

#include <Eigen/Sparse>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <utility>
#include <set>

using namespace mrock::symbolic_operators;
using namespace mrock::symbolic_operators::experimental;
using SparseMatrix = Eigen::SparseMatrix<double>;

constexpr int L = 4;
constexpr int matrix_size = 1 << (2 * L);

SparseMatrix identity(matrix_size,matrix_size);

int wrap_momentum(int momentum) {
    momentum %= L;
    return momentum < 0 ? momentum + L : momentum;
}

double cosine(int momentum) {
    return std::cos(std::numbers::pi * (2. * momentum / L - 1.));
}

SparseMatrix operator_string(const std::vector<int>& modes, const std::vector<bool>& creation) {
    SparseMatrix result(matrix_size, matrix_size);
    if (modes.size() != creation.size()) {
        throw std::invalid_argument("operator and mode counts differ");
    }
    if (modes.empty()) {
        result.setIdentity();
        return result;
    }

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

// We assume that creation is ordered: all true then all false
SparseMatrix normal_ordered_operator_string(const std::vector<int>& modes, const std::vector<bool>& creation, const std::array<double,L>& occupations)
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
        return normal_ordered_operator_string(r_modes, r_creation, occupations);
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

                    //std::cout << i << "," << j << " | " << remaining[0] << "," << remaining[k] 
                    //    << "   -> " << sign1 << " * " << sign2 << " = " << sign1*sign2 << std::endl;

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

TermCollector get_symbolic_H() {
    const Operator c_K_sigma(Momentum('K'), Index::Sigma, false);
    const Operator c_K_sigma_prime(Momentum('K'), Index::SigmaPrime, false);
    return {
        Term(1, Coefficient::RealInversionSymmetric("\\tilde{\\varepsilon}", MomentumList(Momentum('K'))),
             SumContainer{MomentumSum{'K'}, IndexSum{Index::Sigma}},
             {c_K_sigma.hermitian_conjugate(), c_K_sigma}),
        Term(1, Coefficient::RealInteraction("U", MomentumList({Momentum('K'), Momentum('P'), Momentum('Q')}),
                                            IndexWrapper{Index::Sigma, Index::SigmaPrime}),
             SumContainer{MomentumSum{'K', 'P', 'Q'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
             {c_K_sigma.hermitian_conjugate(), c_K_sigma_prime.with_momentum('P').hermitian_conjugate(),
              c_K_sigma_prime.with_momentum(Momentum("P-Q")), c_K_sigma.with_momentum(Momentum("K+Q"))})
    };
}

TermCollector get_symbolic_eta() {
    const Operator c_K_sigma(Momentum('K'), Index::Sigma, false);
    const Operator c_K_sigma_prime(Momentum('K'), Index::SigmaPrime, false);
    return {
        Term(2, Coefficient::RealInteraction("\\alpha", MomentumList({Momentum('K'), Momentum('P'), Momentum('Q')}),
                                                IndexWrapper{Index::Sigma, Index::SigmaPrime}),
                SumContainer{MomentumSum{'K', 'P', 'Q'}, IndexSum{Index::Sigma, Index::SigmaPrime}},
                {
                    c_K_sigma.hermitian_conjugate(), c_K_sigma_prime.with_momentum('P').hermitian_conjugate(),
                    c_K_sigma_prime.with_momentum(Momentum("P-Q")), c_K_sigma.with_momentum(Momentum("K+Q"))
                })
    };
}

std::array<double, L> get_occupation_numbers() {
    std::array<double, L> occupations{};
    constexpr double zero_tolerance = 1.e-12;
    for (int k = 0; k < L; ++k) {
        const double cosine = std::cos(std::numbers::pi * (2. * k / L - 1.));
        occupations[k] = cosine < -zero_tolerance ? 1.
                         : (std::abs(cosine) <= zero_tolerance ? .5 : 0.);
    }
    return occupations;
}

bool matrices_equal(const SparseMatrix& lhs, const SparseMatrix& rhs) {
    return (lhs - rhs).norm() <= 1.e-12;
}

SparseMatrix explicit_total_number_operator() {
    SparseMatrix result(matrix_size, matrix_size);
    for (int momentum = 0; momentum < L; ++momentum) {
        for (int spin = 0; spin < 2; ++spin) {
            const int mode = 2 * momentum + spin;
            result += operator_string({mode, mode}, {true, false});
        }
    }
    return result;
}

SparseMatrix explicit_total_energy() {
    SparseMatrix result(matrix_size, matrix_size);
    for (int momentum = 0; momentum < L; ++momentum) {
        for (int spin = 0; spin < 2; ++spin) {
            const int mode = 2 * momentum + spin;
            result += cosine(momentum) * operator_string({mode, mode}, {true, false});
        }
    }
    return result;
}

template <class MomentumLookup, class IndexLookup>
double coefficient_value(const Coefficient& coefficient, MomentumLookup&& momentum_lookup, IndexLookup&& index_lookup) {
    const int k = momentum_lookup(coefficient.momenta[0]);
    if (coefficient.name == "\\tilde{\\varepsilon}") {
        return cosine(k);
    }
    const int p = momentum_lookup(coefficient.momenta[1]);
    const int q = momentum_lookup(coefficient.momenta[2]);
    const double interaction = cosine(wrap_momentum(k + q)) * cosine(wrap_momentum(p - q));
    if (coefficient.name == "U") {
        return interaction * (index_lookup(coefficient.indices[0]) == index_lookup(coefficient.indices[1]) ? cosine(q) : 1.);
    }
    if (coefficient.name == "\\alpha") {
        return interaction * (cosine(k) + cosine(p) - cosine(wrap_momentum(k + q)) - cosine(wrap_momentum(p - q)))
               * (index_lookup(coefficient.indices[0]) == index_lookup(coefficient.indices[1]) ? cosine(q) : 1.);
    }
    throw std::runtime_error("unknown coefficient: " + coefficient.name);
}

template <class CollectorType>
SparseMatrix symbolic_to_matrix(const CollectorType& terms, [[maybe_unused]] const std::array<double,L>& occupations)
{
    SparseMatrix result(matrix_size, matrix_size);
    for (const auto& term : terms) {
        std::array<int, 256> momenta{};
        std::array<int, 256> indices{};
        const auto momentum_lookup = [&momenta](const Momentum& momentum) {
            int value = 0;
            for (const auto& symbol : momentum.momentum_list) {
                value += symbol.factor * momenta[static_cast<unsigned char>(symbol.name)];
            }
            return wrap_momentum(value + (momentum.add_PI ? L / 2 : 0));
        };
        const auto index_lookup = [&indices](const Index& index) {
            return indices[static_cast<unsigned char>(index)];
        };

        std::function<void(std::size_t)> assign_spins;

        const auto evaluate = [&]() {
            double prefactor = static_cast<double>(term.multiplicity);
            for (const auto& coefficient : term.coefficients) {
                prefactor *= coefficient_value(coefficient, momentum_lookup, index_lookup);
            }

            if constexpr (std::is_same_v<CollectorType, WickOrderedCollector>) {
                for (const auto& expec : term.operators) {
                    if (expec.type != OperatorType::Number) {
                        throw std::runtime_error("Expected only Number-type expectation values!");
                    }
                    prefactor *= occupations[momentum_lookup(expec.momentum)];
                }
            }

            std::vector<int> modes;
            std::vector<bool> creation;

            const std::vector<Operator>& symbolic_operator_string = ([&]() -> const std::vector<Operator>& { 
                if constexpr (std::is_same_v<CollectorType, WickOrderedCollector>) {
                    return term.wick_expression.operators;
                }
                else {
                    return term.operators;
                }
            } )();
            for (const auto& op : symbolic_operator_string) {
                const int spin = is_mutable(op.first_index()) ? indices[static_cast<unsigned char>(op.first_index())]
                                                               : (op.first_index() == Index::SpinDown ? 1 : 0);
                modes.push_back(2 * momentum_lookup(op.momentum) + spin);
                creation.push_back(op.is_daggered);
            }
            if (prefactor != 0.) {
                if constexpr (std::is_same_v<CollectorType, WickOrderedCollector>) {
                    result += prefactor * normal_ordered_operator_string(modes, creation, occupations);
                }
                else {
                    result += prefactor * operator_string(modes, creation);
                }
            }
        };

        assign_spins = [&](std::size_t n) {
            if (n == term.sums.spins.size()) {
                evaluate();
                return;
            }
            indices[static_cast<unsigned char>(term.sums.spins[n])] = 0;
            assign_spins(n + 1);
            indices[static_cast<unsigned char>(term.sums.spins[n])] = 1;
            assign_spins(n + 1);
        };

        std::function<void(std::size_t)> assign_momenta;
        assign_momenta = [&](std::size_t n) {
            if (n == term.sums.momenta.size()) {
                assign_spins(0);
                return;
            }
            const auto name = static_cast<unsigned char>(term.sums.momenta[n]);
            for (int value = 0; value < L; ++value) {
                momenta[name] = value;
                assign_momenta(n + 1);
            }
        };

        assign_momenta(0);
    }
    return result;
}

int main() {
    const std::array<double, L> occupations = get_occupation_numbers();
    const auto templates = NickelCUT::commute::get_wick_templates();
    const auto symmetries = NickelCUT::commute::get_symmetries();
    identity.setIdentity();

    // Verify the algebra
    {
        const SparseMatrix creation     = operator_string({L / 2}, {true});
        const SparseMatrix annihilation = operator_string({L / 2}, {false});
        const SparseMatrix annihilation_transpose = annihilation.transpose();

        double diff = (creation - annihilation_transpose).norm();
        if (diff != 0.) {
            std::cerr << "c^+ != c     diff=" << diff << std::endl;
            return 1;
        }
        diff = (creation * annihilation - operator_string({L / 2, L / 2}, {true, false}))
                   .norm();
        if (diff != 0.) {
            std::cerr << "n != c^+ c     diff=" << diff << std::endl;
            return 1;
        }

        diff = (creation * creation).norm();
        if (diff != 0.) {
            std::cerr << "c c != 0     diff=" << diff << std::endl;
            return 1;
        }
        diff = (creation * annihilation + annihilation * creation).norm() - (matrix_size >> L);
        if (diff != 0.) {
            std::cerr << "{c, c^+} != 1     diff=" << diff << std::endl;
            return 1;
        }

        const SparseMatrix annihilation2 =  operator_string({0}, {false});
        diff = (creation * annihilation2 + annihilation2 * creation).norm();
        if (diff != 0.) {
            std::cerr << "{c_i, c_j^+} != 0     diff=" << diff << std::endl;
            return 1;
        }
        diff = (annihilation * annihilation2 + annihilation2 * annihilation).norm();
        if (diff != 0.) {
            std::cerr << "{c_i^+, c_j^+} != 0     diff=" << diff << std::endl;
            return 1;
        }

        diff = operator_string({0,0,0,0}, {true,true,false,false}).norm();
        if (diff != 0.) {
            std::cerr << "c^+ c^+ c c != 0     diff=" << diff << std::endl;
            return 1;
        }
    }

    // Bilinear tests 1
    {
        for (int K=0; K<L;++K) {
            SparseMatrix n = operator_string({2*K,2*K}, {true,false});
            SparseMatrix no = normal_ordered_operator_string({2*K,2*K}, {true,false}, occupations);
            SparseMatrix expected = n;
            expected -= occupations[K] * identity;

            if ((no-expected).norm() != 0.) {
                std::cerr << "bilinear normal_ordered_operator_string does not yield the expected result!" 
                    << "K=" << K << ": " << (no-expected).norm() << std::endl;
                return 1;
            }
        }
    }

    // Bilinear tests 2
    {
        const Operator c_K_sigma(Momentum('K'), Index::Sigma, false);
        const Term total_number(1, SumContainer{MomentumSum{'K'}, IndexSum{Index::Sigma}},
                                {c_K_sigma.hermitian_conjugate(), c_K_sigma});
        WickOrderedCollector total_number_wick = wick_decompose(total_number, templates);
        clean_wick_ordered_terms(total_number_wick);
        if (!matrices_equal(explicit_total_number_operator(), symbolic_to_matrix(total_number_wick, occupations))) {
            std::cerr << "sum_(k,sigma) c_(k,sigma)^dagger c_(k,sigma) changed under Wick decomposition\n";
            std::cerr << "\\begin{align*}\n" << total_number << "\\end{align*}\n\\begin{align*}" << total_number_wick << "\\end{align*}" << std::endl;  
            return 1;
        }

        const Term total_energy(1, Coefficient("\\tilde{\\varepsilon}", Momentum('K')), 
                                SumContainer{MomentumSum{'K'}, IndexSum{Index::Sigma}},
                                {c_K_sigma.hermitian_conjugate(), c_K_sigma});
        WickOrderedCollector total_energy_wick = wick_decompose(total_energy, templates);
        clean_wick_ordered_terms(total_energy_wick);
        if (!matrices_equal(explicit_total_energy(), symbolic_to_matrix(total_energy_wick, occupations))) {
            std::cerr << "sum_(k,sigma) eps(k) c_(k,sigma)^dagger c_(k,sigma) changed under Wick decomposition\n";
            std::cerr << "\\begin{align*}\n" << total_energy << "\\end{align*}\n\\begin{align*}" << total_energy_wick << "\\end{align*}" << std::endl;  
            return 1;
        }

        WickOrderedCollector total_number_backwards = { WickOrderedTerm() };
        total_number_backwards[0].multiplicity = 1;
        total_number_backwards[0].wick_expression.operators = { c_K_sigma.hermitian_conjugate(), c_K_sigma };
        total_number_backwards[0].sums = SumContainer{MomentumSum{'K'}, IndexSum{Index::Sigma}};
        SparseMatrix backwards_should_be = explicit_total_number_operator();
        double constant_to_substract=0.;
        for (int K=0; K<L;++K) {
            constant_to_substract += 2 * occupations[K];
        }
        backwards_should_be -= constant_to_substract * identity;
        if (!matrices_equal(backwards_should_be, symbolic_to_matrix(total_number_backwards, occupations))) {
            std::cerr << "\\begin{align*}\n" << total_number_backwards << "\\end{align*}" << std::endl;
            std::cerr << "Matrices for total_number_backwards do not match!" << std::endl;
            return 1;
        }
    }
    
    // Quartic tests
    {
        {
            SparseMatrix no = normal_ordered_operator_string({0,0,0,0}, {true,true,false,false}, occupations);
            if ((no).norm() != 0.) {
                std::cerr << "quartic normal_ordered_operator_string does not yield the expected result!" 
                    << "i=" << 0 << ": " << (no).norm() << std::endl;
                return 1;
            }
        }
        for (int i=1; i<2*L;++i) {
            SparseMatrix raw = operator_string({0,i,i,0}, {true,true,false,false})
                            - occupations[0]   * operator_string({i,i},{true, false})
                            - occupations[i/2] * operator_string({0,0},{true, false})
                            + occupations[0]   * occupations[i/2] * identity;
            SparseMatrix no = normal_ordered_operator_string({0,i,i,0}, {true,true,false,false}, occupations);

            if ((no-raw).norm() != 0.) {
                std::cerr << "quartic normal_ordered_operator_string does not yield the expected result!" 
                    << "i=" << i << ": " << (no-raw).norm() << std::endl;
                return 1;
            }
        }
        for (int i=0; i<2*L;++i) {
            if (i == L) continue;
            SparseMatrix raw = operator_string({L,i,i,L}, {true,true,false,false})
                            - occupations[L/2] * operator_string({i,i},{true, false})
                            - occupations[i/2] * operator_string({L,L},{true, false})
                            + occupations[L/2] * occupations[i/2] * identity;
            SparseMatrix no = normal_ordered_operator_string({L,i,i,L}, {true,true,false,false}, occupations);

            if ((no-raw).norm() != 0.) {
                std::cerr << "quartic normal_ordered_operator_string does not yield the expected result!" 
                    << "base=L | i=" << i << ": " << (no-raw).norm() << std::endl;
                return 1;
            }
        }
    }

    const auto check = [&](const TermCollector& expression, const char* name) {
        const SparseMatrix before = symbolic_to_matrix(expression, occupations);
        WickOrderedCollector after = wick_decompose(expression, templates);
        clean_wick_ordered_terms(after, symmetries);
        
        const SparseMatrix after_matrix = symbolic_to_matrix(after, occupations);
        if ((before - after_matrix).norm() > 1.e-12) {
            std::cerr << name << " changed under Wick decomposition. Diff = " << (before - after_matrix).norm() << "\n";
            std::cerr << "\\begin{align*}\n" << expression << "\\end{align*}\n\\begin{align*}" << after << "\\end{align*}" << std::endl;  
            return false;
        }
        return true;
    };

    const TermCollector symbolic_H = get_symbolic_H();
    const TermCollector symbolic_eta = get_symbolic_eta();

    if (!check(symbolic_H, "H")) {
        return 1;
    }
    if (!check(symbolic_eta, "eta")) {
        return 1;
    }

    // sextic tests
    {
        // Degenerate case: all operators on the same mode.
        {
            SparseMatrix no = normal_ordered_operator_string(
                {0,0,0,0,0,0},
                {true,true,true,false,false,false},
                occupations);
            
            if (no.norm() != 0.) {
                std::cerr << "sextic normal_ordered_operator_string does not yield the expected result!"
                          << " i=j=0: " << no.norm() << std::endl;
                return 1;
            }
        }

        // Generic case: three distinct modes 0,i,j.
        for (int i = 1; i < 2*L; ++i) {
            for (int j = i+1; j < 2*L; ++j) {
            
                SparseMatrix raw = operator_string({0,i,j,j,i,0}, {true,true,true,false,false,false})
                    
                    - occupations[0]   * operator_string({i,j,j,i}, {true,true,false,false})
                    - occupations[i/2] * operator_string({0,j,j,0}, {true,true,false,false})
                    - occupations[j/2] * operator_string({0,i,i,0}, {true,true,false,false})
                        
                    + occupations[0]   * occupations[i/2] * operator_string({j,j}, {true,false})
                    + occupations[0]   * occupations[j/2] * operator_string({i,i}, {true,false})
                    + occupations[i/2] * occupations[j/2] * operator_string({0,0}, {true,false})
                        
                    - occupations[0] * occupations[i/2] * occupations[j/2] * identity;
                        
                SparseMatrix no = normal_ordered_operator_string({0,i,j,j,i,0}, {true,true,true,false,false,false}, occupations);
                
                if ((no - raw).norm() != 0.) {
                    std::cerr << "sextic normal_ordered_operator_string does not yield the expected result!"
                              << " i=" << i
                              << " j=" << j
                              << ": " << (no - raw).norm()
                              << std::endl;
                    return 1;
                }
            }
        }
    }

    TermCollector symbolic_commutator = commutator(symbolic_eta, symbolic_H);
    symbolic_commutator.clean_up();

    if (!check(symbolic_commutator, "[eta, H]")) {
        return 1;
    }

    return 0;
}
