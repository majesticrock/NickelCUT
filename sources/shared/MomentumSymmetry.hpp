#pragma once

#include "../L.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <vector>

namespace NickelCUT {

template <std::size_t MomentumCount, std::size_t TupleSize>
class MomentumSymmetry {
    static constexpr std::size_t total_size = [] {
        std::size_t result = 1;
        for (std::size_t i = 0; i < TupleSize; ++i) {
            result *= MomentumCount;
        }
        return result;
    }();
    static constexpr std::size_t fixed_axis_reflection_size = 2 * L;
    static constexpr std::size_t fixed_diagonal_reflection_size = L;
    static constexpr std::size_t fixed_quarter_rotation_size = 2;
    static constexpr std::size_t fixed_half_rotation_size = 4;
    static constexpr std::size_t power(std::size_t base, std::size_t exponent) noexcept {
        std::size_t result = 1;
        for (std::size_t i = 0; i < exponent; ++i) {
            result *= base;
        }
        return result;
    }
    static constexpr std::size_t stored_size =
        (total_size + 2 * power(fixed_quarter_rotation_size, TupleSize) +
         power(fixed_half_rotation_size, TupleSize) + 2 * power(fixed_axis_reflection_size, TupleSize) +
         2 * power(fixed_diagonal_reflection_size, TupleSize)) /
        8;

    static_assert(MomentumCount == static_cast<std::size_t>(N));
    static_assert(stored_size <= static_cast<std::size_t>(std::numeric_limits<int>::max()));

    static std::size_t reflect_kx(std::size_t momentum) noexcept {
        const std::size_t kx = momentum % L;
        const std::size_t ky = momentum / L;
        return (static_cast<std::size_t>(L) - kx) % L + L * ky;
    }

    static std::size_t rotate_90(std::size_t momentum) noexcept {
        const std::size_t kx = momentum % L;
        const std::size_t ky = momentum / L;
        return (static_cast<std::size_t>(L) - ky) % L + L * kx;
    }

    static std::size_t canonical_index(std::size_t index) noexcept {
        std::array<std::size_t, TupleSize> momenta{};
        std::size_t remaining = index;
        for (std::size_t i = TupleSize; i > 0; --i) {
            momenta[i - 1] = remaining % MomentumCount;
            remaining /= MomentumCount;
        }

        std::size_t canonical = index;
        for (int rotation = 0; rotation < 4; ++rotation) {
            std::size_t flattened = 0;
            std::size_t reflected = 0;
            for (std::size_t i = 0; i < TupleSize; ++i) {
                flattened = flattened * MomentumCount + momenta[i];
                reflected = reflected * MomentumCount + reflect_kx(momenta[i]);
            }
            canonical = std::min({canonical, flattened, reflected});
            for (auto& momentum : momenta) {
                momentum = rotate_90(momentum);
            }
        }
        return canonical;
    }

public:
    static constexpr std::size_t independent_size() noexcept { return stored_size; }

    static const std::vector<int>& storage_indices() noexcept {
        static const std::vector<int> indices = [] {
            std::vector<int> result(total_size);
            std::vector<int> canonical_to_storage(total_size, std::numeric_limits<int>::max());

            std::size_t next_storage_index = 0;
            for (std::size_t index = 0; index < total_size; ++index) {
                const std::size_t canonical = canonical_index(index);
                if (canonical == index) {
                    canonical_to_storage[canonical] = static_cast<int>(next_storage_index++);
                }
            }
            assert(next_storage_index == stored_size);

            for (std::size_t index = 0; index < total_size; ++index) {
                result[index] = canonical_to_storage[canonical_index(index)];
            }
            return result;
        }();
        return indices;
    }

    static const std::vector<std::size_t>& independent_indices() noexcept {
        static const std::vector<std::size_t> indices = [] {
            std::vector<std::size_t> result(stored_size);
            std::size_t next_index = 0;
            for (std::size_t index = 0; index < total_size; ++index) {
                if (canonical_index(index) == index) {
                    result[next_index++] = index;
                }
            }
            assert(next_index == stored_size);
            return result;
        }();
        return indices;
    }

    template <class... Indices>
    static int data_index(Indices... indices) noexcept {
        static_assert(sizeof...(Indices) == TupleSize);
        std::size_t flattened = 0;
        ((flattened = flattened * MomentumCount + static_cast<std::size_t>(indices)), ...);
        assert(flattened < total_size);
        return storage_indices()[flattened];
    }
};

}  // namespace NickelCUT
