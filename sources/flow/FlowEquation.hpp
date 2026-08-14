#pragma once

#include "../L.hpp"
#include "FlowContainer.hpp"

#include <array>

namespace NickelCUT::flow {

struct FlowEquation {
	void operator() (const FlowContainer& current, FlowContainer& dHdl, const double /*l*/);
};
} // namespace NickelCUT::flow