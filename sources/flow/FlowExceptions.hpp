#pragma once

#include <stdexcept>
#include <string>

namespace NickelCUT::flow {

class ControlledFlowInterruption : public std::runtime_error {
    const double _end_time;

public:
    ControlledFlowInterruption(const std::string& msg, double end_time);
    ControlledFlowInterruption(const char* msg, double end_time);

    double get_end_time() const noexcept;
};

class LargeRODException : public ControlledFlowInterruption {
public:
    explicit LargeRODException(double end_time);
};

class LongRuntimeException : public ControlledFlowInterruption {
public:
    explicit LongRuntimeException(double end_time);
};

} // namespace NickelCUT::flow
