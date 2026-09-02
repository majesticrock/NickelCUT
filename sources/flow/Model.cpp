#include "Model.hpp"
#include "../helper_functions.hpp"
#include "../L.hpp"

#include <mrock/utility/better_to_string.hpp>
#include <mrock/utility/OutputConvenience.hpp>
#include <nlohmann/json.hpp>

#include <cmath>
#include <string>
#include <sstream>

namespace NickelCUT::flow
{
Model::Model(double U_0_, double tprime_, double E_F_, double temperature_) 
    : U_0{U_0_},
    tprime{tprime_},
    E_F{E_F_},
    temperature{temperature_},
    beta{temperature > 0.0 ? 1. / temperature : -1.}
    //chemical_potential{E_F_}
{}

double Model::epsilon_0(double kx, double ky) const noexcept
{
    return -2. * (
                std::cos(kx) + std::cos(ky) 
                + 2. * tprime * std::cos(kx) * std::cos(ky) 
            ) - E_F;
}

double Model::fermi_function_zero_temperature(double energy) const noexcept 
{
    if (is_zero(energy)) {
        return 0.5;
    }
    return energy < 0.0 ? 1.0 : 0.0;
}

double Model::fermi_function(double energy) const noexcept
{
    if (this->beta >= 0.0) {
        return 1. / (1. + std::exp(beta*energy));
    }
    return fermi_function_zero_temperature(energy);
}

std::string Model::info() const
{
    return "NickelCUT Model (L=" + std::to_string(L) + "): U_0=" + std::to_string(U_0)
        + "   tprime=" + std::to_string(tprime)
        + "   E_F=" + std::to_string(E_F)
        + "   temperature=" + std::to_string(temperature); 
}

std::string Model::data_dir_name() const
{
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

    return "L=" + std::to_string(L) + "/"
        + "T=" + improved_string(temperature) + "/"
        + "U_0=" + improved_string(U_0) + "/"
        + "tprime=" + improved_string(tprime) + "/"
        + "E_F=" + improved_string(E_F) + "/"; 
}

nlohmann::json Model::generate_meta_data_json() const {
    return nlohmann::json{
        { "time", mrock::utility::time_stamp() },
        { "L", L },
        { "T", temperature },
        { "U_0", U_0 },
        { "tprime", tprime },
        { "E_F", E_F },
        { "beta", beta },
        { "filling", filling }
    };
}

} // namespace NickelCUT::flow
