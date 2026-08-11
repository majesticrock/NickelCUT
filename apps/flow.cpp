#include "../sources/flow/L.hpp"
#include "../sources/flow/momentum_iterator.hpp"
#include "../sources/flow/FlowContainer.hpp"
#include "../sources/flow/Model.hpp"
#include "../sources/flow/FlowEquation.hpp"

#include <boost/numeric/odeint.hpp>
#include <mrock/utility/OutputConvenience.hpp>
#include <nlohmann/json.hpp>

#include <list>
#include <iostream>

using namespace NickelCUT::flow;
using namespace boost::numeric::odeint;

int main() {
    Model model{0.5, 0., 0.0, 10};
    FlowContainer flow_state(model);

    FlowEquation flow_equation;

    std::list<double> l_times = {0.};
    std::list<FlowContainer> flow_states = {flow_state};

    FlowContainer dHdl;

    flow_equation(flow_state, dHdl, double{});
    flow_state += 0.01 * dHdl;
    flow_states.push_back(flow_state);

    dHdl = FlowContainer();
    flow_equation(flow_state, dHdl, double{});
    flow_state += 0.01 * dHdl;
    flow_states.push_back(flow_state);

    //constexpr double abs_error = 1e-10;
    //constexpr double rel_error = 1e-6;
    //constexpr double l_final = 0.1;
    //constexpr double dl_output = 0.1;
    //constexpr double dl_use = dl_output / 1;
//
    //double l = 0.;
    //while(l < l_final) {
    //    //integrate_adaptive(make_controlled< runge_kutta_cash_karp54< FlowContainer > >( abs_error, rel_error ),
    //    //            flow_equation, flow_state, l, l + dl_output, dl_use);
//
    //    std::cout << "l=" << l << std::endl;
    //    runge_kutta4<FlowContainer> stepper;
    //    integrate_const(stepper, flow_equation, flow_state, l , l + dl_output , dl_use );
//
    //    l += dl_output;
    //    l_times.push_back(l);
    //    flow_states.push_back(flow_state);
    //}

    nlohmann::json jData = {
        { "time", mrock::utility::time_stamp() },
        { "l_times", l_times },
        { "flow_states", flow_states }
    };

    mrock::utility::save_string_raw(jData.dump(4), "build/test.json");

    return 0;
}