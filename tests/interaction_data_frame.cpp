#include "../sources/flow/InteractionDataFrame.hpp"

#include <iostream>
#include <vector>

using namespace NickelCUT;
using namespace NickelCUT::flow;

namespace {
std::size_t reflect_kx(std::size_t momentum) {
    return (static_cast<std::size_t>(L) - momentum % L) % L + L * (momentum / L);
}

std::size_t reflect_ky(std::size_t momentum) {
    return momentum % L + L * ((static_cast<std::size_t>(L) - momentum / L) % L);
}

std::size_t full_index(std::size_t k, std::size_t p, std::size_t q) {
    return q + N * (p + N * k);
}

std::vector<double> dense_snapshot(const InteractionDataFrame& frame) {
    std::vector<double> result(N * N * N);
    for (std::size_t k = 0; k < N; ++k) {
        for (std::size_t p = 0; p < N; ++p) {
            for (std::size_t q = 0; q < N; ++q) {
                result[full_index(k, p, q)] = frame(k, p, q);
            }
        }
    }
    return result;
}

void fill_distinct(InteractionDataFrame& frame) {
    for (std::size_t index = 0; index < frame.get_data().size(); ++index) {
        frame.get_data()[index] = static_cast<double>(index + 1);
    }
}
}

int main() {
    InteractionDataFrame frame;
    constexpr std::size_t full_size = N * N * N;
    constexpr std::size_t fixed_axis_reflection_size = 2 * L;
    constexpr std::size_t expected_stored_size =
        (full_size + 2 * fixed_axis_reflection_size * fixed_axis_reflection_size * fixed_axis_reflection_size
            + 4 * 4 * 4) / 4;

    if (frame.size() != expected_stored_size) {
        std::cerr << "Unexpected compact size: " << frame.size()
                  << " != " << expected_stored_size << '\n';
        return 1;
    }

    fill_distinct(frame);

    std::vector<bool> seen(frame.size(), false);
    double expected_norm_squared = 0.0;
    for (std::size_t k = 0; k < N; ++k) {
        for (std::size_t p = 0; p < N; ++p) {
            for (std::size_t q = 0; q < N; ++q) {
                const double value = frame(k, p, q);
                const std::size_t stored_index = static_cast<std::size_t>(value) - 1;
                seen[stored_index] = true;
                expected_norm_squared += value * value;

                if (value != frame(reflect_kx(k), reflect_kx(p), reflect_kx(q))
                    || value != frame(reflect_ky(k), reflect_ky(p), reflect_ky(q))) {
                    std::cerr << "A reflected tuple did not share its stored value\n";
                    return 1;
                }
            }
        }
    }

    for (const bool was_seen : seen) {
        if (!was_seen) {
            std::cerr << "A compact storage slot has no logical tuple\n";
            return 1;
        }
    }

    if (frame.abs_squared_total() != expected_norm_squared) {
        std::cerr << "The norm did not account for reflected tuple multiplicity\n";
        return 1;
    }

    frame(1, 3, 5) = -7.0;
    if (frame(reflect_kx(1), reflect_kx(3), reflect_kx(5)) != -7.0
        || frame(reflect_ky(1), reflect_ky(3), reflect_ky(5)) != -7.0) {
        std::cerr << "Writing through a reflected tuple did not update its representative\n";
        return 1;
    }

    InteractionDataFrame symmetrized;
    fill_distinct(symmetrized);
    const auto before_symmetrize = dense_snapshot(symmetrized);
    symmetrized.symmetrize();
    for (std::size_t k = 0; k < N; ++k) {
        for (std::size_t p = 0; p < N; ++p) {
            for (std::size_t q = 0; q < N; ++q) {
                const momentum_iterator<L> K(static_cast<int>(k));
                const momentum_iterator<L> P(static_cast<int>(p));
                const momentum_iterator<L> Q(static_cast<int>(q));
                const double expected = 0.5 * (
                    before_symmetrize[full_index(k, p, q)]
                    + before_symmetrize[full_index(P, K, -Q)]
                );
                if (symmetrized(k, p, q) != expected) {
                    std::cerr << "symmetrize produced an unexpected value\n";
                    return 1;
                }
            }
        }
    }

    InteractionDataFrame antisymmetrized;
    fill_distinct(antisymmetrized);
    const auto before_antisymmetrize = dense_snapshot(antisymmetrized);
    antisymmetrized.antisymmetrize();
    for (std::size_t k = 0; k < N; ++k) {
        for (std::size_t p = 0; p < N; ++p) {
            for (std::size_t q = 0; q < N; ++q) {
                const momentum_iterator<L> K(static_cast<int>(k));
                const momentum_iterator<L> P(static_cast<int>(p));
                const momentum_iterator<L> Q(static_cast<int>(q));
                const double expected = 0.5 * (
                    before_antisymmetrize[full_index(k, p, q)]
                    - before_antisymmetrize[full_index(K, P, P - K - Q)]
                );
                if (antisymmetrized(k, p, q) != expected) {
                    std::cerr << "antisymmetrize produced an unexpected value\n";
                    return 1;
                }
            }
        }
    }

    return 0;
}