#pragma once

#include "FlowContainer.hpp"

namespace NickelCUT::flow {

struct FlowEquation {
	void operator() (const FlowContainer& current, FlowContainer& dHdl, const double /*l*/);
};
} // namespace NickelCUT::flow