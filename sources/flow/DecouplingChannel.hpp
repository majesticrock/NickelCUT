#pragma once

#include "FlowContainer.hpp"

#include <cassert>
#include <cstddef>
#include <vector>
#include <utility>

namespace NickelCUT::flow
{

class DecouplingChannel {  
    const int _N; // _L*_L
    const int total_size; // _N*_N
    
    // Storage order: q is contiguous because q is the sum variable in the mean-field equations
    std::vector<double> _data;

public:
    template<class Archive>
    void serialize(Archive& ar, [[maybe_unused]] const unsigned int version) {
        ar & _data;
    }

    DecouplingChannel() = default;
    DecouplingChannel(int N);

    // The first element is for differing spins, the second for parallel spins
    static std::pair<DecouplingChannel, DecouplingChannel> SingleParticleEnergy(const FlowContainer& base);
    // The first element is for differing spins, the second for parallel spins
    static std::pair<DecouplingChannel, DecouplingChannel> DensityWave(const FlowContainer& base);
    // Only differing spins contribute
    static DecouplingChannel Superconductivity(const FlowContainer& base); 

    inline double& operator()(std::size_t p, std::size_t q) noexcept {
        assert(p < static_cast<std::size_t>(_N));
        assert(q < static_cast<std::size_t>(_N));
        return _data[_N*p + q];
    }
    inline double operator()(std::size_t p, std::size_t q) const noexcept {
        assert(p < static_cast<std::size_t>(_N));
        assert(q < static_cast<std::size_t>(_N));
        return _data[_N*p + q];
    }

    std::vector<std::vector<double>> as_2D_array() const noexcept;
}; 

} // namespace NickelCUT::flow
