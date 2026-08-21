#include "DecouplingChannel.hpp"

#include "../L.hpp"

namespace NickelCUT::flow
{

DecouplingChannel::DecouplingChannel(int N) 
    : _N(N), total_size(_N*_N), _data(total_size*total_size)
{}

std::pair<DecouplingChannel, DecouplingChannel> DecouplingChannel::SingleParticleEnergy(const FlowContainer& base)
{
    using mom_it = FlowContainer::mom_it;

    std::pair<DecouplingChannel, DecouplingChannel> ret = { DecouplingChannel(N), DecouplingChannel(N) };

    for (mom_it p = mom_it::begin(); p != mom_it::end(); ++p) {
        for (mom_it q = mom_it::begin(); q != mom_it::end(); ++q) {
            ret.first(p, q)  = 2. * base.interactions_differing_spin(p, q, Gamma<L>);
            ret.second(p, q) = 2. * (
                base.interactions_same_spin(p, q, Gamma<L>) - base.interactions_same_spin(p, q, q - p)
            );
        }
    }
    return ret;
}

std::pair<DecouplingChannel, DecouplingChannel> DecouplingChannel::DensityWave(const FlowContainer& base)
{
    using mom_it = FlowContainer::mom_it;

    std::pair<DecouplingChannel, DecouplingChannel> ret = { DecouplingChannel(N), DecouplingChannel(N) };

    for (mom_it p = mom_it::begin(); p != mom_it::end(); ++p) {
        for (mom_it q = mom_it::begin(); q != mom_it::end(); ++q) {
            ret.first(p, q)  = 2. * base.interactions_differing_spin(p, q, Q<L>);
            ret.second(p, q) = 2. * (
                base.interactions_same_spin(p, q, Q<L>) - base.interactions_same_spin(p, q, q - p + Q<L>)
            );
        }
    }
    return ret;
}

DecouplingChannel DecouplingChannel::Superconductivity(const FlowContainer& base)
{
    using mom_it = FlowContainer::mom_it;

    DecouplingChannel ret(N);

    for (mom_it p = mom_it::begin(); p != mom_it::end(); ++p) {
        for (mom_it q = mom_it::begin(); q != mom_it::end(); ++q) {
            ret(p, q) = 2. * base.interactions_differing_spin(p, -p, q - p);
        }
    }
    return ret;
}

std::vector<std::vector<double>> DecouplingChannel::as_2D_array() const noexcept {
    std::vector<std::vector<double>> ret(N, std::vector<double>(N));
    for (auto p = momentum_iterator<L>::begin(); p != momentum_iterator<L>::end(); ++p) {
        for (auto q = momentum_iterator<L>::begin(); q != momentum_iterator<L>::end(); ++q) {
            ret[p][q] = operator()(p, q);
        }
    }
    return ret;
}

} // namespace NickelCUT::flow
