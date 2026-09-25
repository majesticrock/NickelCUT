#pragma once

#include "../L.hpp"
#include "FlowContainer.hpp"
#include <boost/numeric/odeint.hpp>

namespace NickelCUT::flow
{
inline constexpr double _ROD_0(double U_0) {
    return 0.5 * (U_0 < 0. ? -1. : 1.) * U_0;
}

inline constexpr double l_final(double U_0) {
    return (U_0 < 0. ? -1. : 1.) * (10. / U_0);
}

inline constexpr double target_dl(double U_0) {
    return l_final(U_0) / 50.;
}

inline constexpr double dl(double U_0) {
    return target_dl(U_0) / 50.;
}

constexpr double abs_error = 1e-6;
constexpr double rel_error = 1e-6;

typedef boost::numeric::odeint::runge_kutta_fehlberg78<FlowContainer, double, FlowContainer, double, boost::numeric::odeint::vector_space_algebra> boost_stepper;
} // namespace NickelCUT::flow

// Required for boost::odeint controlled steppers
namespace boost { namespace numeric { namespace odeint {
template<>
struct vector_space_norm_inf< NickelCUT::flow::FlowContainer >
{
    typedef double result_type;
    double operator()( const NickelCUT::flow::FlowContainer &p ) const
    {
        return p.norm_inf();
    }
};
} } } // boost namespaces