#include "TCFinder.hpp"
#include "../helper_functions.hpp"

#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>

namespace NickelCUT::mean_field
{

constexpr double MEAN_FIELD_PRECISION = 1e-10;

void to_json(nlohmann::json& j, const TransitionData& transition) noexcept
{
    j = nlohmann::json{
        { "lower_temperature", transition.lower_temperature },
        { "upper_temperature", transition.upper_temperature },
        { "lower_filling",     transition.lower_filling     },
        { "upper_filling",     transition.upper_filling     },
        { "order",             to_string(transition.order)  }
    };
    j.update(transition.sc_result);
}

TCFinder::TCFinder(Model& _model, double _dT) 
    : model(_model), 
    initial_sc_values(_model.deltas),
    dT(_dT),
    sc_solver(&model, &model.deltas, 200, 1e-6)
{}

void TCFinder::bracket_phase_transistion(TransitionData& transition)
{
    constexpr double temperature_precision = 1e-4;

    while (transition.upper_temperature - transition.lower_temperature > temperature_precision) {
        model.temperature = 0.5 * (transition.lower_temperature + transition.upper_temperature);
        model.beta = 1. / model.temperature;
        reset_model();

        sc_solver.compute(false, 600, MEAN_FIELD_PRECISION);
        if (!model.deltas.converged) break;

        if (model.order_type(MEAN_FIELD_PRECISION) == transition.order) {
            transition.lower_temperature = model.temperature;
            transition.lower_filling = model.compute_filling();
        } else {
            transition.upper_temperature = model.temperature;
            transition.upper_filling = model.compute_filling();
        }
    }

    model.temperature = transition.upper_temperature;
    model.beta = 1. / model.temperature;
}

void TCFinder::reset_model() 
{
    model.deltas = initial_sc_values;
    const double max = *std::max_element(initial_sc_values.begin(), initial_sc_values.begin() + 3*N, LessThanAbs());
    for (int i=0; i<3*N; ++i) {
        model.deltas[i] += 0.1 * max * std::cos(0.12345678*i);
    }
}

std::deque<TransitionData> TCFinder::compute()
{
    std::deque<TransitionData> ret;
    OrderType current_order = model.order_type(MEAN_FIELD_PRECISION);
    OrderType new_order;
    double old_temperature;
    double old_filling;

    while (current_order != OrderType::NormalMetal) {
        old_temperature = model.temperature;
        old_filling = model.compute_filling();

        model.temperature += dT;
        model.beta = 1. / model.temperature;

        reset_model();
        sc_solver.compute(false, 600, MEAN_FIELD_PRECISION);
        if (!model.deltas.converged) {
            std::cout << "T_c iteration did not converge... Retrying..." << std::endl;
            sc_solver.compute(false, 600, MEAN_FIELD_PRECISION);
            if (!model.deltas.converged) {
                std::cerr << "T_c iteration did not converge after retrying..." << std::endl;
                break;
            }
        }

        new_order = model.order_type(MEAN_FIELD_PRECISION);

        if (new_order != current_order) {
            // A phase transistion has occurred!
            TransitionData transition{
                old_temperature, 
                model.temperature, 
                old_filling, 
                model.compute_filling(), 
                current_order,
                model.selfconsistency_to_json()
            };
            bracket_phase_transistion(transition);
            ret.push_back(transition);
        }

        current_order = new_order;

        std::cout << "Computed model for T=" << model.temperature << " and found " << to_string(current_order) << std::endl;
    }

    return ret;
}

} // namespace NickelCUT::mean_field
