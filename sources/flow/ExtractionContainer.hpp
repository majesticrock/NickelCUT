#pragma once

#include "FlowContainer.hpp"
#include "DecouplingChannel.hpp"

#include <nlohmann/json.hpp>

#include <utility>

namespace NickelCUT::flow {

struct ExtractionContainer {
    // The first element is for differing spins, the second for parallel spins
    std::pair<DecouplingChannel, DecouplingChannel> single_particle_energy;
    // The first element is for differing spins, the second for parallel spins
    std::pair<DecouplingChannel, DecouplingChannel> density_wave;
    // Only differing spins contribute
    DecouplingChannel superconductivity;

    FlowContainer::coeff_array dispersion;

    FlowContainer::coeff_array epsilon_tilde;

    template<class Archive>
    void serialize(Archive& ar, [[maybe_unused]] const unsigned int version) {
        ar & single_particle_energy;
        ar & density_wave;
        ar & superconductivity;
        ar & dispersion;
        ar & epsilon_tilde;
    }

    ExtractionContainer(const FlowContainer& x);
    ExtractionContainer() noexcept = default;
};

void to_json(nlohmann::json& j, const ExtractionContainer& extracted_channels) noexcept;

} // namespace NickelCUT::flow