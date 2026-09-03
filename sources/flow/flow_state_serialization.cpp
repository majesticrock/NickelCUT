#include "flow_state_serialization.hpp"
#include "FlowContainer.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>

#include <string>
#include <fstream>
#include <stdexcept>

namespace NickelCUT::flow
{

void serialize_flow_state(const FlowContainer& state, const std::string& output_dir, const std::string& output_filename)
{
    const std::string file = output_dir + output_filename;
    std::ofstream ofs(file, std::ios::binary);
    if (ofs.good()) {
        boost::archive::binary_oarchive oa(ofs);
        oa << state;
    }
    else {
        throw std::runtime_error("Outputstream for " + file + " is bad!");
    }
}

FlowContainer deserialize_flow_state(const std::string& input_dir, const std::string& input_filename)
{
    const std::string file = input_dir + input_filename;
    std::ifstream ifs(file, std::ios::binary);
    if (ifs.good()) {
        FlowContainer state;

        boost::archive::binary_iarchive ia(ifs);
        ia >> state;

        return state;
    }
    else {
        throw std::runtime_error("Inputstream for " + file + " is bad!");
    }
}

} // namespace NickelCUT::flow
