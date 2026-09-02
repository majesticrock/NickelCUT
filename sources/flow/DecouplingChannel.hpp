#pragma once

#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

namespace NickelCUT::flow
{

struct FlowContainer;

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

    DecouplingChannel() : _N(0), total_size(0) {}
    DecouplingChannel(int N);

    DecouplingChannel(const DecouplingChannel&) = default;
    DecouplingChannel(DecouplingChannel&& other) noexcept
        : _N(other._N), total_size(other.total_size), _data(std::move(other._data)) {}

    DecouplingChannel& operator=(const DecouplingChannel& other) {
        assert(_N == other._N);
        assert(total_size == other.total_size);
        _data = other._data;
        return *this;
    }
    DecouplingChannel& operator=(DecouplingChannel&& other) noexcept {
        assert(_N == other._N);
        assert(total_size == other.total_size);
        _data = std::move(other._data);
        return *this;
    }

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
