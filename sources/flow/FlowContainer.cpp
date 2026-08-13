#include "FlowContainer.hpp"
#include "momentum_iterator.hpp"
#include "compare_abs.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>

namespace NickelCUT::flow {

FlowContainer::FlowContainer() 
    : interactions_same_spin(), interactions_differing_spin()
{
    dispersion.fill(double{});
    occupation_numbers.fill(double{});
    epsilon_tilde.fill(double{});
}

FlowContainer::FlowContainer(const Model& model) 
    : interactions_same_spin(), interactions_differing_spin(model.U_0 / N)
{
    for (mom_it p = mom_it::begin(); p != mom_it::end(); ++p) {
        dispersion[p.get_position()] = model.epsilon_0(p.get_kx(), p.get_ky());
        occupation_numbers[p.get_position()] = model.fermi_function(dispersion[p.get_position()]);
    }
    fill_epsilon_tilde();
}

void FlowContainer::fill(double value) noexcept
{
    interactions_same_spin.get_data().fill(value);
    interactions_differing_spin.get_data().fill(value);
    dispersion.fill(value);
    occupation_numbers.fill(value);
    epsilon_tilde.fill(value);
}

void FlowContainer::reset() noexcept
{
    fill(double{});
}

bool FlowContainer::contains_nan_or_inf() const noexcept {
    for (const auto& val : interactions_same_spin.get_data()) {
        if (std::isnan(val) || std::isinf(val)) return true;
    }
    for (const auto& val : interactions_differing_spin.get_data()) {
        if (std::isnan(val) || std::isinf(val)) return true;
    }
    for (const auto& val : dispersion) {
        if (std::isnan(val) || std::isinf(val)) return true;
    }
    for (const auto& val : epsilon_tilde) {
        if (std::isnan(val) || std::isinf(val)) return true;
    }
    return false;
}

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
        { "dispersion",                  container.dispersion                                },
        //{ "epsilon_tilde",               container.epsilon_tilde                             },
        { "interactions_same_spin",      container.interactions_same_spin.as_3D_array()      },
        { "interactions_differing_spin", container.interactions_differing_spin.as_3D_array() }
    };
}

} // namespace NickelCUT::flow