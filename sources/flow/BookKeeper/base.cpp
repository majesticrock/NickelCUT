#include "base.hpp"

#include "../ExtractionContainer.hpp"
#include "../FlowEquation.hpp"
#include "../FlowExceptions.hpp"
#include "../../shared/FlowContainer.hpp"

#include <mrock/utility/OutputConvenience.hpp>

#include <deque>
#include <iostream>

namespace NickelCUT::flow::BookKeeper {
base::base(const FlowContainer& initial_flow_state,
           FlowEquation& flow_equation,
           double _band_width,
           double _dl,
           std::chrono::minutes::rep _max_runtime_duration)
    : l_times{0.0},
      residual_offdiagonalities{initial_flow_state.residual_offdiagonality()},
      max_interactions{N * initial_flow_state.max_interaction()},
      last_good_state{initial_flow_state},
      flow_equation{flow_equation},
      dl{_dl},
      initial_band_width{_band_width},
      max_runtime_duration{_max_runtime_duration},
      begin(clock::now()),
      last(begin),
      current_idx{0U} {
    derivative_residual_offdiagonalities.push_back(compute_derivative_rod(initial_flow_state, 0.0));
    std::cout << mrock::utility::time_stamp() << "   -   " << "Starting calculations...\n"
              << "Initial ROD = " << residual_offdiagonalities.back() << "\n"
              << "Initial derivative ROD = " << derivative_residual_offdiagonalities.back() << "\n"
              << "Saving data with a spacing of at least dl=" << dl << std::endl;
};

double base::compute_derivative_rod(const FlowContainer& state, double l) const {
    FlowContainer derivative;
    flow_equation(state, derivative, l);
    return derivative.residual_offdiagonality();
}

nlohmann::json base::convert_to_json() const noexcept {
    return nlohmann::json{{"l_times", l_times},
                          {"number_of_data_points", l_times.size()},
                          {"residual_offdiagonalities", residual_offdiagonalities},
                          {"derivative_residual_offdiagonalities", derivative_residual_offdiagonalities},
                          {"max_interactions", max_interactions}};
}

void base::basic_processing(const FlowContainer& x, const double l) {
    const auto total_runtime = std::chrono::duration_cast<std::chrono::minutes>(clock::now() - begin).count();
    if (total_runtime > max_runtime_duration) {
        throw LongRuntimeException(l);
    }
    if (l - l_times.back() < dl)
        return;

    const bool state_is_good = x.is_inversion_symmetric() && x.is_hermitian();
    if (!state_is_good) {
        throw BrokenSymmetriesException(l_times.back());
    }

    const double max_coeff = N * std::max(x.interactions_differing_spin.norm_inf(), x.interactions_same_spin.norm_inf());
    const double ROD = x.residual_offdiagonality();
    const double derivative_ROD = compute_derivative_rod(x, l);

    if (max_coeff > 10. * initial_band_width) {
        throw LargeInteractionException(l);
    }
    if (derivative_ROD < 0.002 * initial_band_width) {
        throw SmallDerivativeRODException(l);
    }
    //if (derivative_ROD > 5. * initial_band_width) {
    //    throw LargeDerivativeRODException(l);
    //}
    
    ++current_idx;
    clock::time_point now = clock::now();
    std::cout << "//------------------------------------------------------//\n"
              << "Step #" << current_idx << "\t" << mrock::utility::time_stamp() << "\n"
              << "Step took " << std::chrono::duration_cast<std::chrono::milliseconds>(now - last).count()
              << "ms to execute." 
              << "l = " << l << "\t\tmax(U) = " << max_coeff 
              << "\t\tROD = " << ROD 
              << "\t\tROD(dH/dl) = " << derivative_ROD << "\n"
              << std::endl;
    last = now;

    l_times.push_back(l);
    residual_offdiagonalities.push_back(ROD);
    derivative_residual_offdiagonalities.push_back(derivative_ROD);
    max_interactions.push_back(max_coeff);
    last_good_state = x;
}

bool base::basic_wrap_up(const FlowContainer& x, const double l) {
    const bool accept_last_step = x.is_inversion_symmetric() && x.is_hermitian()
        && !float_equal(l_times.back(), l);

    if (accept_last_step) {
        const double derivative_rod = compute_derivative_rod(x, l);
        l_times.push_back(l);
        residual_offdiagonalities.push_back(x.residual_offdiagonality());
        derivative_residual_offdiagonalities.push_back(derivative_rod);
        max_interactions.push_back(N * x.max_interaction());
        last_good_state = x;
    }

    const std::chrono::hh_mm_ss hms(clock::now() - begin);

    std::cout << "//------------------------------------------------------//\n"
              << "\t Flow program finished at " << mrock::utility::time_stamp() << "\n"
              << "\t\tfinal max(U) = " << max_interactions.back() << "\n"
              << "Total executation took " << std::setfill('0') << std::setw(2) << hms.hours().count() << ':'
              << std::setw(2) << hms.minutes().count() << ':' << std::setw(2) << hms.seconds().count() << '\n'
              << "\tGoodbye." << std::endl;

    return accept_last_step;
}

}  // namespace NickelCUT::flow::BookKeeper
