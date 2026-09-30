#include "ExtractionContainer.hpp"
#include "DecouplingChannel.hpp"

#include <nlohmann/json.hpp>

#include <stdexcept>
#include <string>
#include <vector>

namespace NickelCUT::flow 
{

ExtractionContainer::ExtractionContainer(const FlowContainer& x)
    : single_particle_energy(DecouplingChannel::SingleParticleEnergy(x)),
    density_wave(DecouplingChannel::DensityWave(x)),
    superconductivity(DecouplingChannel::Superconductivity(x)),
    dispersion(x.dispersion),
    epsilon_tilde(x.epsilon_tilde)
{ }

void to_json(nlohmann::json& j, const ExtractionContainer& extracted_channels) noexcept
{
    j = nlohmann::json{
        { "single_particle_energy_differing",   extracted_channels.single_particle_energy.first.as_2D_array()  },
        { "single_particle_energy_same",        extracted_channels.single_particle_energy.second.as_2D_array() },
        { "density_wave_differing",             extracted_channels.density_wave.first.as_2D_array()            },
        { "density_wave_same",                  extracted_channels.density_wave.second.as_2D_array()           },
        { "superconductivity",                  extracted_channels.superconductivity.as_2D_array()             },
        { "epsilon_tilde",                      extracted_channels.epsilon_tilde                               },
        { "dispersion",                         extracted_channels.dispersion                                  }
    };
}

void from_json(const nlohmann::json& j, ExtractionContainer& extracted_channels)
{
    const auto read_channel = [&j](const char* key) {
        const auto matrix = j.at(key).get<std::vector<std::vector<double>>>();
        if (matrix.size() != N) {
            throw std::invalid_argument(std::string("Invalid channel row count for ") + key);
        }

        DecouplingChannel channel(N);
        for (std::size_t row = 0; row < N; ++row) {
            if (matrix[row].size() != N) {
                throw std::invalid_argument(std::string("Invalid channel column count for ") + key);
            }
            for (std::size_t column = 0; column < N; ++column) {
                channel(row, column) = matrix[row][column];
            }
        }
        return channel;
    };

    extracted_channels.single_particle_energy.first = read_channel("single_particle_energy_differing");
    extracted_channels.single_particle_energy.second = read_channel("single_particle_energy_same");
    extracted_channels.density_wave.first = read_channel("density_wave_differing");
    extracted_channels.density_wave.second = read_channel("density_wave_same");
    extracted_channels.superconductivity = read_channel("superconductivity");
    extracted_channels.dispersion = j.at("dispersion").get<FlowContainer::coeff_array>();
    extracted_channels.epsilon_tilde = j.at("epsilon_tilde").get<FlowContainer::coeff_array>();
}

} // namespace NickelCUT::flow