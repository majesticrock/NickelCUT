#pragma once

#include "../L.hpp"
#include "FlowContainer.hpp"
#include <boost/numeric/odeint.hpp>

namespace NickelCUT::flow
{
constexpr double U_0 = -1.;
constexpr double T = 0.;
constexpr double tprime = 0.;
constexpr double mu_0 = 0.01;

constexpr double _ROD_0 = 0.5 * (U_0 < 0. ? -1. : 1.) * U_0 * L;

constexpr double l_final   = (50. / _ROD_0);
constexpr double target_dl = (1. / (5. * L * _ROD_0));
constexpr double dl = target_dl / 50.;

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