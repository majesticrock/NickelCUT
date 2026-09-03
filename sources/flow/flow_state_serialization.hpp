#pragma once

#include <string>

namespace NickelCUT::flow
{

struct FlowContainer;

void serialize_flow_state(const FlowContainer& state, const std::string& output_dir, const std::string& output_filename);

FlowContainer deserialize_flow_state(const std::string& input_dir, const std::string& input_filename);

} // namespace NickelCUT::flow
