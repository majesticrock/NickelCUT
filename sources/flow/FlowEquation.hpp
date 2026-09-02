#pragma once

namespace NickelCUT::flow {

struct FlowContainer;

struct FlowEquation {
	void operator() (const FlowContainer& current, FlowContainer& dHdl, const double /*l*/);
};
} // namespace NickelCUT::flow