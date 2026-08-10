#pragma once

#include <mrock/symbolic_operators/WickTerm.hpp>

namespace NickelCUT::commute {
using namespace mrock::symbolic_operators;

void export_as_flow_equation(const WickTermCollector& bilinears, const WickTermCollector& quartics);

} // namespace NickelCUT::commute