#pragma once

#include "momentum_iterator.hpp"

#include <array>
#include <algorithm>
#include <iostream>

namespace NickelCUT::flow {

namespace detail {
    template<int _L>
    constexpr static int _N_uniques = (_L/2 + 1)*(_L/2 + 2)/2;

    template<int _L>
    constexpr std::array<momentum_iterator<_L>, _N_uniques<_L>> make_uniques() noexcept {
        std::array<momentum_iterator<_L>, _N_uniques<_L>> ret;
        int i=0;
        for (int y=0; y<=_L/2; ++y) {
            for (int x=y; x<=_L/2; ++x) {
                ret[i] = momentum_iterator<_L>(x, y);
                ++i;
            }
        }
        return ret;
    }
}

template<int _L>
struct unique_momenta {
    static_assert(_L % 2 == 0);
    constexpr static int N_uniques = detail::_N_uniques<_L>;

    constexpr static std::array<momentum_iterator<_L>, N_uniques> uniques = detail::make_uniques<_L>();

    constexpr static int reduce(momentum_iterator<_L> momentum) noexcept {
        if (momentum._y == _L/2) {
            if (momentum._x == _L/2) return momentum;
            if (momentum._x > _L/2) momentum._x = _L - momentum._x;
            momentum._y = momentum._x;
            momentum._x = _L/2;
            return momentum_iterator<_L>(momentum._x, momentum._y);
        }
        if (momentum._y > _L/2) {
            momentum = -momentum;
        }

        if (momentum._x > _L/2) {
            momentum._x = _L - momentum._x;
        }
        // y and x are both < _L/2
        if (momentum._x < momentum._y) {
            std::swap(momentum._x, momentum._y);
        }

        return momentum_iterator<_L>(momentum._x, momentum._y);
    }

    template<class UnaryContainerAccess>
    constexpr static void promote(UnaryContainerAccess&& container) noexcept(UnaryContainerAccess()) {
        for (momentum_iterator<_L> k = momentum_iterator<_L>::begin(); k != momentum_iterator<_L>::end(); ++k) {
            container(k) = container(reduce(k));
        }
    }
};

} // namespace NickelCUT::flow