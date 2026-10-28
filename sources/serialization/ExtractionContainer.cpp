#include "ExtractionContainer.hpp"

#include "../flow/ExtractionContainer.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <fstream>
#include <stdexcept>
#include <string>

namespace NickelCUT::flow {

void serialize_extracted_channels(const ExtractionContainer& state,
                                  const std::string& output_dir,
                                  const std::string& output_filename) {
    const std::string file = output_dir + output_filename;
    std::ofstream ofs(file, std::ios::binary);
    if (ofs.good()) {
        boost::archive::binary_oarchive oa(ofs);
        oa << state;
    } else {
        throw std::runtime_error("Outputstream for " + file + " is bad!");
    }
}

ExtractionContainer deserialize_extracted_channels(const std::string& input_dir, const std::string& input_filename) {
    const std::string file = input_dir + input_filename;
    std::ifstream ifs(file, std::ios::binary);
    if (ifs.good()) {
        ExtractionContainer extracted;

        boost::archive::binary_iarchive ia(ifs);
        ia >> extracted;

        return extracted;
    } else {
        throw std::runtime_error("Inputstream for " + file + " is bad!");
    }
}

}  // namespace NickelCUT::flow
