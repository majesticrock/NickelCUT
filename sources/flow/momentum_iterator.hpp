#pragma once

#include <cassert>
#include <array>
#include <numbers>
#include <ostream>
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
struct momentum_iterator;

template <int _L>
static constexpr momentum_iterator<_L> Gamma = momentum_iterator<_L>(_L / 2, _L / 2);

template <int _L>
static constexpr momentum_iterator<_L> Q = momentum_iterator<_L>(0, 0);

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
        assert(n >= 0 && n < _L*_L);
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
    
    constexpr momentum_iterator() = default;

    constexpr momentum_iterator(int x, int y) noexcept : _x(x), _y(y), _pos(x + _L*y) {}

    constexpr momentum_iterator(int pos) noexcept : _x(pos % _L), _y(pos / _L), _pos(pos) {}

    constexpr operator std::size_t() const noexcept {
        return static_cast<std::size_t>(_pos);
    }

    static constexpr momentum_iterator begin() noexcept {
        return momentum_iterator(0, 0);
    }
    static constexpr momentum_iterator end() noexcept {
        return momentum_iterator(0, _L);
    }
    // For iterating half of the first BZ
    static constexpr momentum_iterator half_end() noexcept {
        return momentum_iterator(0, _L / 2);
    }

    constexpr bool operator==(const momentum_iterator& other) const noexcept { return _pos == other._pos; }
    constexpr bool operator!=(const momentum_iterator& other) const noexcept { return !(*this == other); }

    constexpr momentum_iterator& operator+=(momentum_iterator other) noexcept {
        _x += other._x + 3 * _L / 2;
        _y += other._y + 3 * _L / 2;
        assert(_x >= 0);
        assert(_y >= 0);
        _x %= _L;
        _y %= _L;
        _pos = _x + _L*_y;
        return *this;
    }
    constexpr momentum_iterator& operator-=(momentum_iterator other) noexcept {
        _x -= other._x - 3 * _L / 2;
        _y -= other._y - 3 * _L / 2; 
        assert(_x >= 0);
        assert(_y >= 0);
        _x %= _L;
        _y %= _L;
        _pos = _x + _L*_y;
        return *this;
    }
    constexpr momentum_iterator& operator*=(int factor) noexcept {
        if (factor == 0) {
            *this = Gamma<_L>;
            return *this;
        }

        const momentum_iterator original = *this;
        *this = Gamma<_L>;
        if (factor > 0) {
            for (int count = 0; count < factor; ++count) {
                *this += original;
            }
        } else {
            for (int count = 0; count > factor; --count) {
                *this -= original;
            }
        }
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
    constexpr momentum_iterator operator*(int factor) const noexcept {
        auto tmp = *this;
        tmp *= factor;
        return tmp;
    }
    constexpr momentum_iterator operator-() const noexcept {
        return Gamma<_L> - (*this);
    }
    
private:
    int _x{};
    int _y{};
    int _pos{};
};

template<int _L>
std::ostream& operator<<(std::ostream& os, momentum_iterator<_L> mom_it) {
    os << "(" << mom_it.get_kx() / std::numbers::pi << ", " << mom_it.get_ky() / std::numbers::pi << ")";
    return os; 
}

template<int _L>
constexpr momentum_iterator<_L> operator*(int factor, momentum_iterator<_L> momentum) noexcept {
    momentum *= factor;
    return momentum;
}

static_assert(Gamma<6>.get_kx() == 0.0);
static_assert(Gamma<6>.get_ky() == 0.0);

static_assert(Gamma<10>.get_kx() == 0.0);
static_assert(Gamma<10>.get_ky() == 0.0);

} // namespace NickelCUT::flow