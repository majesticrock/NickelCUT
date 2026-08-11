#pragma once

#include <array>
#include <numbers>
#include <utility>

namespace NickelCUT::flow {

namespace detail {
    template<int _L>
    constexpr std::array<double, _L> fill_momentum_cache() {
        std::array<double, _L> ret;
        for (int i=0; i<_L; ++i) {
            ret[i] = std::numbers::pi * (static_cast<double>(2 * i) / static_cast<double>(_L) - 1.0);
        }
        return ret;
    }
}

template<int _L>
struct momentum_iterator {
    static_assert(_L % 2 == 0);
    static constexpr std::array<double, _L> momentum_cache = detail::fill_momentum_cache<_L>();

    constexpr double get_kx() const noexcept {
        return momentum_cache[_x];
    }
    constexpr double get_ky() const noexcept {
        return momentum_cache[_y];
    }

    constexpr int get_position() const noexcept {
        return _pos;
    }

    constexpr std::pair<double, double> get_momentum() const noexcept {
        return {momentum_cache[_x], momentum_cache[_y]};
    }

    constexpr std::pair<double, double> operator[](int n) const noexcept {
        return {momentum_cache[n % _L], momentum_cache[n / _L]};
    }

    constexpr momentum_iterator& operator++() noexcept {
        if(++_x >= _L) {
            _x = 0;
            ++_y;
        } 
        ++_pos;
        return *this;
    }

    constexpr momentum_iterator operator++(int) noexcept {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    constexpr momentum_iterator(int x, int y) noexcept : _x(x), _y(y), _pos(x + _L*y) {}

    static constexpr momentum_iterator begin() noexcept {
        return momentum_iterator(0, 0);
    }
    static constexpr momentum_iterator end() noexcept {
        return momentum_iterator(0, _L);
    }

    constexpr bool operator==(const momentum_iterator& other) const noexcept { return _pos == other._pos; }
    constexpr bool operator!=(const momentum_iterator& other) const noexcept { return !(*this == other); }

    constexpr momentum_iterator& operator+=(momentum_iterator other) noexcept {
        _x += other._x;
        _y += other._y;
        _x %= _L;
        _y %= _L;
        _pos = _x + _L*_y;
        return *this;
    }
    constexpr momentum_iterator& operator-=(momentum_iterator other) noexcept {
        _x -= other._x - _L; // + _L to avoid negative values; _L will be cancelled by the modulo anyway
        _y -= other._y - _L; 
        _x %= _L;
        _y %= _L;
        _pos = _x + _L*_y;
        return *this;
    }
    
    constexpr momentum_iterator operator+(momentum_iterator other) const noexcept {
        auto tmp = *this;
        tmp += other;
        return tmp;
    }
    constexpr momentum_iterator operator-(momentum_iterator other) const noexcept {
        auto tmp = *this;
        tmp -= other;
        return tmp;
    }
    constexpr momentum_iterator& operator-() noexcept {
        _x = _L - _x;
        _y = _L - _y;
        return *this;
    }
private:
    int _x{};
    int _y{};
    int _pos{};
};

template <int _L>
static constexpr momentum_iterator<_L> GammaPoint = momentum_iterator<_L>(_L / 2, _L / 2);

static_assert(GammaPoint<10>.get_kx() == 0.0);
static_assert(GammaPoint<10>.get_ky() == 0.0);

} // namespace NickelCUT::flow