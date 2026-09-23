#pragma once

#include "Model.hpp"
#include "OrderType.hpp"

#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>
#include <nlohmann/json.hpp>

#include <deque>

namespace NickelCUT::mean_field
{

struct TransitionData {
    double lower_temperature; ///< The highest temperature before the phase transition
    double upper_temperature; ///< The lowest temperature after the phase transistion
    // This means, the phase transistion occurs somewhere between lower_temperature and upper_temperature
    double lower_filling; ///< Filling of the system at the lower temperature
    double upper_filling; ///< Filling of the system at the upper temperature
    // The filling cannot properly be kept constant due to the coarse grid
    // Thus, we save both fillings to also have a bracket in the filling
    OrderType order; ///< The order before the phase transistion (i.e. at low temperatures)

    nlohmann::json sc_result;
};

void to_json(nlohmann::json& j, const TransitionData& transition) noexcept;


class TCFinder {
    Model& model;
    decltype(model.deltas) initial_sc_values;
    double dT;

    typedef mrock::utility::Selfconsistency::BroydenSolver<double, Model, decltype(model.deltas), mrock::utility::Selfconsistency::NoWarning> sc_solver_t;
    sc_solver_t sc_solver;

    void bracket_phase_transistion(TransitionData& transition);

    void reset_model();
public:
    TCFinder(Model& _model, double _dT);

    // Assumes that the mean-field calculations on model have run before
    std::deque<TransitionData> compute();
};

} // namespace NickelCUT::mean_field
