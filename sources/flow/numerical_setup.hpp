#pragma once

#include "../L.hpp"
#include "../shared/FlowContainer.hpp"
#include "../shared/Model.hpp"

#include <boost/numeric/odeint.hpp>

namespace NickelCUT::flow {
inline constexpr double _ROD_0(double U_0) {
    return 0.5 * (U_0 < 0. ? -1. : 1.) * U_0;
}

double finite_size_energy_gap(const Model& model);

inline constexpr double target_dl(double l_final) {
    return l_final / 25.;
}

inline constexpr double dl(double U_0) {
    return target_dl(U_0) / 50.;
}

constexpr double abs_error = 1e-6;
constexpr double rel_error = 1e-6;

typedef boost::numeric::odeint::
    runge_kutta_fehlberg78<FlowContainer, double, FlowContainer, double, boost::numeric::odeint::vector_space_algebra>
        boost_stepper;
}  // namespace NickelCUT::flow

// Required for boost::odeint controlled steppers
namespace boost {
namespace numeric {
namespace odeint {
template <>
struct vector_space_norm_inf<NickelCUT::FlowContainer> {
    typedef double result_type;
    double operator()(const NickelCUT::FlowContainer& p) const { return p.norm_inf(); }
};
}  // namespace odeint
}  // namespace numeric
}  // namespace boost