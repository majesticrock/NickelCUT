#include "FlowContainer.hpp"
#include "momentum_iterator.hpp"
#include "compare_abs.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>

namespace NickelCUT::flow {

void FlowContainer::fill_epsilon_tilde() {
    for (mom_it p = mom_it::begin(); p != mom_it::end(); ++p) {
        epsilon_tilde[p.get_position()] = 0;

        for (mom_it q = mom_it::begin(); q != mom_it::end(); ++q) {
            epsilon_tilde[p.get_position()] += (interactions_differing_spin(p.get_position(), q.get_position(), GammaPoint<L>.get_position())
                                                   + interactions_same_spin(p.get_position(), q.get_position(), GammaPoint<L>.get_position())
                                                   - interactions_same_spin(p.get_position(), q.get_position(), (p-q).get_position())
                    
                                                ) * occupation_numbers[q.get_position()];
        }
        epsilon_tilde[p.get_position()] *= -2.0;
        epsilon_tilde[p.get_position()] += dispersion[p.get_position()];
    }
};

double FlowContainer::abs() const
{
    double val = interactions_same_spin.abs() + interactions_differing_spin.abs();
    for (const auto& element : dispersion) {
        val += std::abs(element);
    }
    return val;
}
double FlowContainer::norm_inf() const
{
    double val = interactions_same_spin.norm_inf();

    double comp = interactions_differing_spin.norm_inf();
    if(val < comp) val = comp;

    comp = *std::max_element(dispersion.begin(), dispersion.end(), LessThanAbs());
    if(val < comp) val = comp;

    return val;
}

FlowContainer& FlowContainer::operator+=(const FlowContainer& other)
{
    interactions_same_spin += other.interactions_same_spin;
    interactions_differing_spin += other.interactions_differing_spin;
    for (int i=0; i<N; ++i) {
        dispersion[i] += other.dispersion[i];
    }

    fill_epsilon_tilde();

    return *this;
}
FlowContainer& FlowContainer::operator-=(const FlowContainer& other)
{
    interactions_same_spin -= other.interactions_same_spin;
    interactions_differing_spin -= other.interactions_differing_spin;
    for (int i=0; i<N; ++i) {
        dispersion[i] -= other.dispersion[i];
    }

    fill_epsilon_tilde();

    return *this;
}

FlowContainer& FlowContainer::operator*=(const FlowContainer& other)
{
    interactions_same_spin *= other.interactions_same_spin;
    interactions_differing_spin *= other.interactions_differing_spin;
    for (int i=0; i<N; ++i) {
        dispersion[i] *= other.dispersion[i];
    }

    fill_epsilon_tilde();

    return *this;
}
FlowContainer& FlowContainer::operator/=(const FlowContainer& other)
{
    interactions_same_spin /= other.interactions_same_spin;
    interactions_differing_spin /= other.interactions_differing_spin;
    for (int i=0; i<N; ++i) {
        dispersion[i] /= other.dispersion[i];
    }

    fill_epsilon_tilde();

    return *this;
}

FlowContainer& FlowContainer::operator*=(const double other)
{
    interactions_same_spin *= other;
    interactions_differing_spin *= other;
    for (int i=0; i<N; ++i) {
        dispersion[i] *= other;
    }

    fill_epsilon_tilde();

    return *this;
}
FlowContainer& FlowContainer::operator/=(const double other)
{
    interactions_same_spin /= other;
    interactions_differing_spin /= other;
    for (int i=0; i<N; ++i) {
        dispersion[i] /= other;
    }

    fill_epsilon_tilde();

    return *this;
}

void to_json(nlohmann::json& j, const FlowContainer& container)
{
    j = nlohmann::json{
        { "dispersion",                  container.dispersion                             },
        //{ "epsilon_tilde",               container.epsilon_tilde                          },
        { "interactions_same_spin",      container.interactions_same_spin.get_data()      },
        { "interactions_differing_spin", container.interactions_differing_spin.get_data() }
    };
}

} // namespace NickelCUT::flow