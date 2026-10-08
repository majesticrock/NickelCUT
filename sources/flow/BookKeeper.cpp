#include "BookKeeper.hpp"
#include "DecouplingChannel.hpp"
#include "FlowContainer.hpp"
#include "FlowExceptions.hpp"
#include "ExtractionContainer.hpp"
#include "FlowEquation.hpp"

#include <mrock/utility/OutputConvenience.hpp>

#include <utility>
#include <deque>
#include <iostream>

namespace NickelCUT::flow
{
BookKeeper::BookKeeper(const FlowContainer& initial_flow_state, FlowEquation& flow_equation, double _band_width, double _dl, std::chrono::minutes::rep _max_runtime_duration)
    : l_times{ 0.0 },
    residual_offdiagonalities{ initial_flow_state.residual_offdiagonality() },
    max_interactions{ N * initial_flow_state.max_interaction() },
    last_good_state{ initial_flow_state },
    flow_equation{ flow_equation },
    dl{ _dl },
    initial_band_width{ _band_width },
    max_runtime_duration{ _max_runtime_duration },
    begin(clock::now()), 
    last(begin),
    current_idx{ 0U }
{
    derivative_residual_offdiagonalities.push_back(compute_derivative_rod(initial_flow_state, 0.0));
    std::cout << mrock::utility::time_stamp() << "   -   " << "Starting calculations...\n"
        << "Initial ROD = " << residual_offdiagonalities.back() << "\n"
        << "Initial derivative ROD = " << derivative_residual_offdiagonalities.back() << "\n"
        << "Saving data with a spacing of at least dl=" << dl << std::endl;
};

double BookKeeper::compute_derivative_rod(const FlowContainer& state, double l) {
    FlowContainer derivative;
    flow_equation(state, derivative, l);
    return derivative.residual_offdiagonality();
}

void BookKeeper::print_final(const FlowContainer& x, double l) {
    const bool state_is_good = x.is_inversion_symmetric() && x.is_hermitian();

    if (!float_equal(l_times.back(), l) && state_is_good) {
        const double derivative_rod = compute_derivative_rod(x, l);
        l_times.push_back(l);
        residual_offdiagonalities.push_back(x.residual_offdiagonality());
        derivative_residual_offdiagonalities.push_back(derivative_rod);
        max_interactions.push_back(N * x.max_interaction());
        last_good_state = x;
    }

    const std::chrono::hh_mm_ss hms(clock::now() - begin);

    std::cout << "//------------------------------------------------------//\n"
        << "\t Flow program finished at "
        << mrock::utility::time_stamp() << "\n"
        << "\t\tfinal max(U) = " << max_interactions.back() << "\n"
        << "Total executation took " 
            << std::setfill('0')
              << std::setw(2) << hms.hours().count() << ':'
              << std::setw(2) << hms.minutes().count() << ':'
              << std::setw(2) << hms.seconds().count()
              << '\n'
        << "\tGoodbye."
        << std::endl;
}

void BookKeeper::operator()(const FlowContainer &x, double l)
{
    const auto total_runtime = std::chrono::duration_cast<std::chrono::minutes>(clock::now() - begin).count();
    if (total_runtime > max_runtime_duration) {
        throw LongRuntimeException(l);
    }
    const double max_coeff = N * std::max(x.interactions_differing_spin.norm_inf(), x.interactions_same_spin.norm_inf());
    if (max_coeff > 10. * initial_band_width) {
        throw LargeInteractionException(l);
    }
    const bool state_is_good = x.is_inversion_symmetric() && x.is_hermitian();
    if (!state_is_good) {
        throw BrokenSymmetriesException(l_times.back());
    }
    const double derivative_rod = compute_derivative_rod(x, l);
    if (derivative_rod > 5. * initial_band_width) {
        throw LargeDerivativeRODException(l);
    }
    if (derivative_rod < 0.001 * initial_band_width) {
        throw SmallDerivativeRODException(l);
    }
    if (l - l_times.back() < dl) return;

    ++current_idx;
    clock::time_point now = clock::now();
    std::cout << "//------------------------------------------------------//\n"
        << "Step #" << current_idx << "\t" << mrock::utility::time_stamp() << "\n"
        << "l = " << l
        << "\t\tmax(U) = " << max_coeff << "\t\tROD(dH/dl) = " << derivative_rod << "\n"
        << "Step took " << std::chrono::duration_cast<std::chrono::milliseconds>(now - last).count() << "ms to execute."
        << std::endl;
    last = now;

    l_times.push_back(l);
    residual_offdiagonalities.push_back(x.residual_offdiagonality());
    derivative_residual_offdiagonalities.push_back(derivative_rod);
    max_interactions.push_back(max_coeff);
    last_good_state = x;
}

void to_json(nlohmann::json& j, const BookKeeper& book_keeper) noexcept
{
    j = nlohmann::json{
        { "l_times",                    book_keeper.l_times                              },
        { "number_of_data_points",      book_keeper.l_times.size()                       },
        { "residual_offdiagonalities",  book_keeper.residual_offdiagonalities            },
        { "derivative_residual_offdiagonalities", book_keeper.derivative_residual_offdiagonalities },
        { "max_interactions",           book_keeper.max_interactions                     },
        { "extracted_channels",         ExtractionContainer(book_keeper.last_good_state) }
    };
}

} // namespace NickelCUT::flow
