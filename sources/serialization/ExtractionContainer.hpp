#pragma once

#include <string>

namespace NickelCUT::flow {

struct ExtractionContainer;

void serialize_extracted_channels(const ExtractionContainer& state,
                                  const std::string& output_dir,
                                  const std::string& output_filename);
ExtractionContainer deserialize_extracted_channels(const std::string& input_dir, const std::string& input_filename);

}  // namespace NickelCUT::flow
