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

LargeInteractionException::LargeInteractionException(double end_time) 
    : ControlledFlowInterruption("One interaction matrix element grew very large...", end_time) {}

BrokenSymmetriesException::BrokenSymmetriesException(double end_time) 
    : ControlledFlowInterruption("Numerical errors have acrued and symmetries are broken.", end_time) {}

}  // namespace NickelCUT::flow
