#pragma once

#include <string>

namespace NickelCUT::flow
{

struct FlowContainer;
struct ExtractionContainer;

void serialize_flow_state(const FlowContainer& state, const std::string& output_dir, const std::string& output_filename);
FlowContainer deserialize_flow_state(const std::string& input_dir, const std::string& input_filename);

void serialize_extracted_channels(const ExtractionContainer& state, const std::string& output_dir, const std::string& output_filename);
ExtractionContainer deserialize_extracted_channels(const std::string& input_dir, const std::string& input_filename);

} // namespace NickelCUT::flow
