#pragma once
namespace NickelCUT {
struct FlowContainer;

namespace flow {
struct FlowEquation {
    void operator()(const FlowContainer& current, FlowContainer& dHdl, const double /*l*/);
};
}  // namespace flow
}  // namespace NickelCUT