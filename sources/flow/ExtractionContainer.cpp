#include "ExtractionContainer.hpp"
#include "DecouplingChannel.hpp"

#include <nlohmann/json.hpp>

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

} // namespace NickelCUT::flow