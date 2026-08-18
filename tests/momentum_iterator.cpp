#include "../sources/flow/momentum_iterator.hpp"

#include <iostream>
#include <cmath>
#include <numbers>

using namespace NickelCUT::flow;
using std::numbers::pi;

constexpr int L = 6;
using it_type = momentum_iterator<L>;

double momentum_modulo(double input) {
    while(input < -pi) input += 2 * pi;
    // - 1e-12 because floating point precision
    while(input >= pi - 1e-12) input -= 2 * pi;
    return input;
}

bool float_equal(double l, double r) {
    return std::abs(l - r) < 1e-12;
}

int main(int /*argc*/, char** /*argv*/) {
    for(int y = 0; y < L; ++y) {
        for(int x = 0; x < L; ++x) {
            it_type from_coordinates(x, y);
            it_type from_position(x + L * y);

            if(from_coordinates.get_position() != x + L * y
                || from_position.get_position() != x + L * y
                || from_coordinates.get_momentum() != from_position.get_momentum()) {
                std::cerr << "coordinate and position constructors disagree!" << std::endl;
                return 1;
            }

            const auto indexed_momentum = it_type::begin()[x + L * y];
            if(indexed_momentum != from_coordinates.get_momentum()) {
                std::cerr << "operator[] does not use row-major positions!" << std::endl;
                return 1;
            }
        }
    }

    if(Q<L>.get_position() != 0 || Q<L>.get_momentum() != it_type::begin().get_momentum()) {
        std::cerr << "Q is not the first momentum point!" << std::endl;
        return 1;
    }

    if(GammaPoint<L>.get_position() != L / 2 + L * (L / 2)
        || GammaPoint<L>.get_momentum() != std::pair<double, double>{0.0, 0.0}) {
        std::cerr << "GammaPoint is not at zero momentum!" << std::endl;
        return 1;
    }

    int half_count = 0;
    for(it_type p = it_type::begin(); p != it_type::half_end(); ++p) {
        ++half_count;
    }
    if(half_count != L * L / 2 || it_type::half_end().get_position() != L * L / 2) {
        std::cerr << "half_end does not cover half of the Brillouin zone!" << std::endl;
        return 1;
    }

    {
        it_type p(4);
        const it_type before = p++;
        if(before.get_position() != 4 || p.get_position() != 5) {
            std::cerr << "postfix increment has incorrect semantics!" << std::endl;
            return 1;
        }

        const it_type& same = ++p;
        if(&same != &p || p.get_position() != 6) {
            std::cerr << "prefix increment has incorrect semantics!" << std::endl;
            return 1;
        }
    }

    for(it_type p = it_type::begin(); p != it_type::end(); ++p) {
        for(it_type add = it_type::begin(); add != it_type::end(); ++add) {
            it_type in_place_add = p;
            in_place_add += add;
            if(in_place_add != p + add) {
                std::cerr << "operator+= disagrees with operator+!" << std::endl;
                return 1;
            }

            it_type in_place_subtract = p;
            in_place_subtract -= add;
            if(in_place_subtract != p - add) {
                std::cerr << "operator-= disagrees with operator-!" << std::endl;
                return 1;
            }
        }
    }

    int n = 0;
    for(it_type p = it_type::begin(); p != it_type::end(); ++p) {
        if(p.get_position() != n) {
            std::cerr << "p.get_position() != n     failed!" << std::endl;
            return 1;
        }

        it_type tmp = p + GammaPoint<L>;
        if(p != tmp){
            std::cerr << "p + GammaPoint<L> == p     failed!" << std::endl;
            return 1;
        }

        it_type tmp2 = p - GammaPoint<L>;
        if(p != tmp2){
            std::cerr << "p - GammaPoint<L> == p     failed!" << std::endl;
            return 1;
        }

        ++n;
    }

    for(it_type p = it_type::begin(); p != it_type::end(); ++p) {
        for (it_type add = it_type::begin(); add != it_type::end(); ++add) {
            it_type tmp = p + add;
            double should_be_x = momentum_modulo(add.get_kx() + p.get_kx());
            double should_be_y = momentum_modulo(add.get_ky() + p.get_ky());
            
            if(!float_equal(should_be_x, tmp.get_kx())){
                std::cerr << "p_x + q_x == p_x + q_x     failed!" << std::endl;
                std::cerr << "Got " << tmp.get_kx() << ", but expected " 
                    << p.get_kx() << " + " << add.get_kx() << " = " << should_be_x << std::endl;
                std::cerr << "Difference is " << should_be_x - tmp.get_kx() << std::endl;
                return 1;
            }
            if(!float_equal(should_be_y, tmp.get_ky())){
                std::cerr << "p_y + q_y == p_y + q_y     failed!" << std::endl;
                std::cerr << "Got " << tmp.get_ky() << ", but expected " 
                    << p.get_ky() << " + " << add.get_ky() << " = " << should_be_y << std::endl;
                std::cerr << "Difference is " << should_be_y - tmp.get_ky() << std::endl;
                return 1;
            }
        }
    }

    for(it_type p = it_type::begin(); p != it_type::end(); ++p) {
        for (it_type add = it_type::begin(); add != it_type::end(); ++add) {
            it_type tmp = p - add;
            double should_be_x = momentum_modulo(p.get_kx() - add.get_kx());
            double should_be_y = momentum_modulo(p.get_ky() - add.get_ky());

            if(!float_equal(should_be_x, tmp.get_kx())){
                std::cerr << "p_x - q_x == p_x - q_x     failed!" << std::endl;
                std::cerr << "Got " << tmp.get_kx() << ", but expected " 
                    << p.get_kx() << " + " << add.get_kx() << " = " << should_be_x << std::endl;
                std::cerr << "Difference is " << should_be_x - tmp.get_kx() << std::endl;
                return 1;
            }
            if(!float_equal(should_be_y, tmp.get_ky())){
                std::cerr << "p_y - q_y == p_y - q_y     failed!" << std::endl;
                std::cerr << "Got " << tmp.get_ky() << ", but expected " 
                    << p.get_ky() << " - " << add.get_ky() << " = " << should_be_y << std::endl;
                std::cerr << "Difference is " << should_be_y - tmp.get_ky() << std::endl;
                return 1;
            }
        }
    }

    for(it_type p = it_type::begin(); p != it_type::end(); ++p) {
        it_type tmp = -p;
        double should_be_x = momentum_modulo(- (p.get_kx()));
        double should_be_y = momentum_modulo(- (p.get_ky()));

        if(!float_equal(should_be_x, tmp.get_kx())){
            std::cerr << "-p_x = -p_x     failed!" << std::endl;
            std::cerr << "Got " << tmp.get_kx() << ", but expected -" << p.get_kx() << " = " << should_be_x << std::endl;
            std::cerr << "Difference is " << should_be_x - tmp.get_kx() << std::endl;
            return 1;
        }
        if(!float_equal(should_be_y, tmp.get_ky())){
            std::cerr << "-p_y = -p_y     failed!" << std::endl;
            std::cerr << "Got " << tmp.get_ky() << ", but expected -" << p.get_ky() << " = " << should_be_y << std::endl;
            std::cerr << "Difference is " << should_be_y - tmp.get_ky() << std::endl;
            return 1;
        }
    }

    return 0;
}