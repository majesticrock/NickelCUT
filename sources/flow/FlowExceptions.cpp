#include "FlowExceptions.hpp"

namespace NickelCUT::flow {

ControlledFlowInterruption::ControlledFlowInterruption(const std::string& msg, double end_time)
    : std::runtime_error(msg), _end_time(end_time) {}

ControlledFlowInterruption::ControlledFlowInterruption(const char* msg, double end_time)
    : std::runtime_error(msg), _end_time(end_time) {}

double ControlledFlowInterruption::get_end_time() const noexcept {
    return _end_time;
}

LargeRODException::LargeRODException(double end_time)
    : ControlledFlowInterruption("The ROD grew very large...", end_time) {}

LongRuntimeException::LongRuntimeException(double end_time)
    : ControlledFlowInterruption("The program ran for a long time...", end_time) {}

} // namespace NickelCUT::flow
