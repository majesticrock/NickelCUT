#include "../sources/shared/InteractionDataFrame.hpp"
#include "../sources/shared/SingleParticleDataFrame.hpp"

#include <iostream>
#include <vector>

using namespace NickelCUT;

namespace {
std::size_t reflect_kx(std::size_t momentum) {
    return (static_cast<std::size_t>(L) - momentum % L) % L + L * (momentum / L);
}

std::size_t reflect_ky(std::size_t momentum) {
    return momentum % L + L * ((static_cast<std::size_t>(L) - momentum / L) % L);
}

std::size_t rotate_90(std::size_t momentum) {
    const std::size_t kx = momentum % L;
    const std::size_t ky = momentum / L;
    return (static_cast<std::size_t>(L) - ky) % L + L * kx;
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
}  // namespace

int main() {
    InteractionDataFrame frame;
    constexpr std::size_t full_size = N * N * N;
    constexpr std::size_t fixed_axis_reflection_size = 2 * L;
    constexpr std::size_t fixed_diagonal_reflection_size = L;
    constexpr std::size_t fixed_quarter_rotation_size = 2;
    constexpr std::size_t fixed_half_rotation_size = 4;
    constexpr std::size_t expected_stored_size =
        (full_size + 2 * fixed_quarter_rotation_size * fixed_quarter_rotation_size * fixed_quarter_rotation_size +
         fixed_half_rotation_size * fixed_half_rotation_size * fixed_half_rotation_size +
         2 * fixed_axis_reflection_size * fixed_axis_reflection_size * fixed_axis_reflection_size +
         2 * fixed_diagonal_reflection_size * fixed_diagonal_reflection_size * fixed_diagonal_reflection_size) /
        8;

    if (frame.size() != expected_stored_size) {
        std::cerr << "Unexpected compact size: " << frame.size() << " != " << expected_stored_size << '\n';
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

                if (value != frame(reflect_kx(k), reflect_kx(p), reflect_kx(q)) ||
                    value != frame(reflect_ky(k), reflect_ky(p), reflect_ky(q)) ||
                    value != frame(rotate_90(k), rotate_90(p), rotate_90(q))) {
                    std::cerr << "A reflected or rotated tuple did not share its stored value\n";
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
    if (frame(reflect_kx(1), reflect_kx(3), reflect_kx(5)) != -7.0 ||
        frame(reflect_ky(1), reflect_ky(3), reflect_ky(5)) != -7.0 ||
        frame(rotate_90(1), rotate_90(3), rotate_90(5)) != -7.0) {
        std::cerr << "Writing through a reflected or rotated tuple did not update its representative\n";
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
                const double expected =
                    0.5 * (before_symmetrize[full_index(k, p, q)] + before_symmetrize[full_index(P, K, -Q)]);
                if (symmetrized(k, p, q) != expected) {
                    std::cerr << "symmetrize produced an unexpected value\n";
                    return 1;
                }
            }
        }
    }

    SingleParticleDataFrame single_particle;
    constexpr std::size_t single_particle_size = (N + 2 * 2 + 4 + 2 * 2 * L + 2 * L) / 8;
    if (single_particle.size() != single_particle_size) {
        std::cerr << "Unexpected single-particle compact size: " << single_particle.size() << " != "
                  << single_particle_size << '\n';
        return 1;
    }

    for (std::size_t index = 0; index < single_particle.get_data().size(); ++index) {
        single_particle.get_data()[index] = static_cast<double>(index + 1);
    }

    std::vector<bool> single_particle_seen(single_particle.size(), false);
    double expected_single_particle_norm_squared = 0.0;
    for (std::size_t momentum = 0; momentum < N; ++momentum) {
        const double value = single_particle(momentum);
        const std::size_t storage_index = static_cast<std::size_t>(value) - 1;
        single_particle_seen[storage_index] = true;
        expected_single_particle_norm_squared += value * value;

        if (value != single_particle(reflect_kx(momentum)) || value != single_particle(reflect_ky(momentum)) ||
            value != single_particle(rotate_90(momentum))) {
            std::cerr << "A reflected or rotated single-particle momentum did not share its stored value\n";
            return 1;
        }
    }

    for (const bool was_seen : single_particle_seen) {
        if (!was_seen) {
            std::cerr << "A single-particle storage slot has no logical momentum\n";
            return 1;
        }
    }

    if (single_particle.abs_squared_total() != expected_single_particle_norm_squared) {
        std::cerr << "The single-particle norm did not account for symmetry multiplicity\n";
        return 1;
    }

    single_particle(1) = -7.0;
    if (single_particle(reflect_kx(1)) != -7.0 || single_particle(reflect_ky(1)) != -7.0 ||
        single_particle(rotate_90(1)) != -7.0) {
        std::cerr << "Writing through a reflected or rotated momentum did not update its representative\n";
        return 1;
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
                const double expected = 0.5 * (before_antisymmetrize[full_index(k, p, q)] -
                                               before_antisymmetrize[full_index(K, P, P - K - Q)]);
                if (antisymmetrized(k, p, q) != expected) {
                    std::cerr << "antisymmetrize produced an unexpected value\n";
                    return 1;
                }
            }
        }
    }

    return 0;
}