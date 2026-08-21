#pragma once

#include "../../experimental/WickOrderedCollector.hpp"

#include <sstream>

namespace NickelCUT::commute::Hamiltonian
{
    mrock::symbolic_operators::experimental::WickOrderedCollector commute_and_normal_order(std::ostringstream& oss);
    
    void commute_and_extract(std::ostringstream& oss);
} // namespace NickelCUT::commute::Hamiltonian
