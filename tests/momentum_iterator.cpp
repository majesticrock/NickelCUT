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

int main() {
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