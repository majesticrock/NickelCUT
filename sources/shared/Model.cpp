#include "Model.hpp"

#include "../L.hpp"
#include "../helper_functions.hpp"
#include "../momentum_iterator.hpp"
#include "occupation_numbers.hpp"

#include <boost/math/tools/toms748_solve.hpp>
#include <mrock/utility/InputFileReader.hpp>
#include <mrock/utility/OutputConvenience.hpp>
#include <mrock/utility/better_to_string.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace NickelCUT {
Model::Model(mrock::utility::InputFileReader& input)
    : U_0{input.getDouble("U_0")},
      tprime{input.getDouble("tprime")},
      filling{input.getDouble("filling")},
      temperature{input.getDouble("T")},
      beta{temperature > 0.0 ? 1. / temperature : 100.} {
    E_F = compute_chemical_potential();
    compute_occupation_numbers(*this);
    std::cout << info() << std::endl;
}

double Model::compute_chemical_potential() const {
    if (!std::isfinite(filling) || filling <= 0.0 || filling >= 2.0) {
        throw std::invalid_argument("The target filling must be finite and between 0 and 2.");
    }

    double min_energy = std::numeric_limits<double>::infinity();
    double max_energy = -std::numeric_limits<double>::infinity();
    for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
        const double energy = epsilon_0(p.get_kx(), p.get_ky()) + E_F;
        min_energy = std::min(min_energy, energy);
        max_energy = std::max(max_energy, energy);
    }

    const double filling_per_spin = filling / 2.0;
    const double logit = std::log(filling_per_spin) - std::log1p(-filling_per_spin);
    const double thermal_shift = logit / beta;
    const double lower_bound = min_energy + thermal_shift;
    const double upper_bound = max_energy + thermal_shift;

    const auto filling_deviation = [this](double chemical_potential) {
        double computed_filling{};
        for (momentum_iterator<L> p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
            const double energy = epsilon_0(p.get_kx(), p.get_ky()) + E_F - chemical_potential;
            computed_filling += fermi_function(energy);
        }
        computed_filling *= 2.0 / N;
        return filling - computed_filling;
    };

    if (filling_deviation(lower_bound) == 0.0) {
        return lower_bound;
    }
    if (filling_deviation(upper_bound) == 0.0) {
        return upper_bound;
    }

    std::uintmax_t max_iterations{100U};
    const auto bounds = boost::math::tools::toms748_solve(
        filling_deviation, lower_bound, upper_bound, boost::math::tools::eps_tolerance<double>(36), max_iterations);

    return 0.5 * (bounds.first + bounds.second);
}

double Model::epsilon_0(double kx, double ky) const noexcept {
    return -2. * (std::cos(kx) + std::cos(ky) + 2. * tprime * std::cos(kx) * std::cos(ky)) - E_F;
}

double Model::fermi_function_zero_temperature(double energy) const noexcept {
    if (is_zero(energy)) {
        return 0.5;
    }
    return energy < 0.0 ? 1.0 : 0.0;
}

double Model::fermi_function(double energy) const noexcept {
    if (this->beta >= 0.0) {
        return 1. / (1. + std::exp(beta * energy));
    }
    return fermi_function_zero_temperature(energy);
}

double Model::band_width() const noexcept {
    double min = 100000.;
    double max = -100000.;
    for (const auto& kx : momentum_iterator<L>::momentum_cache) {
        for (const auto& ky : momentum_iterator<L>::momentum_cache) {
            const double E = epsilon_0(kx, ky);
            if (E < min)
                min = E;
            if (E > max)
                max = E;
        }
    }
    return max - min;
}

std::string Model::info() const {
    return "NickelCUT Model (L=" + std::to_string(L) + "): U_0=" + std::to_string(U_0) +
           "   tprime=" + std::to_string(tprime) + "   filling=" + std::to_string(filling) +
           "   E_F=" + std::to_string(E_F) +
           "   temperature=" + std::to_string(temperature);
}

std::string Model::data_dir_name() const {
    auto improved_string = [](double number) -> std::string {
        if (std::floor(number) == number) {
            // If the number is a whole number, format it with one decimal place
            std::ostringstream out;
            out.precision(1);
            out << std::fixed << number;
            return out.str();
        } else {
            std::string str = mrock::utility::better_to_string(number, std::chars_format::fixed);
            // Remove trailing zeroes
            str.erase(str.find_last_not_of('0') + 1, std::string::npos);
            str.erase(str.find_last_not_of('.') + 1, std::string::npos);
            return str;
        }
    };

    return "L=" + std::to_string(L) + "/" + "T=" + improved_string(temperature) + "/" + "U_0=" + improved_string(U_0) +
           "/" + "tprime=" + improved_string(tprime) + "/" + "filling=" + improved_string(filling) + "/";
}

nlohmann::json Model::generate_meta_data_json() const {
    return nlohmann::json{{"time", mrock::utility::time_stamp()},
                          {"L", L},
                          {"T", temperature},
                          {"U_0", U_0},
                          {"tprime", tprime},
                          {"E_F", E_F},
                          {"beta", beta},
                          {"filling", filling}};
}

}  // namespace NickelCUT
