#pragma once

#include <mrock/utility/InputFileReader.hpp>
#include <nlohmann/json.hpp>

#include <string>

namespace NickelCUT {

struct Model {
    const double U_0{};
    const double tprime{};
    double filling{};
    double temperature{};

    double beta{};
    double E_F{};

    Model(mrock::utility::InputFileReader& input);

    double epsilon_0(double kx, double ky) const noexcept;

    double fermi_function_zero_temperature(double energy) const noexcept;

    double fermi_function(double energy) const noexcept;

    double band_width() const noexcept;

    std::string info() const;

    std::string data_dir_name() const;

    nlohmann::json generate_meta_data_json() const;

private:
    double compute_chemical_potential() const;
};

}  // namespace NickelCUT
