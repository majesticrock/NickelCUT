#include "../sources/commute/Hamiltonian/commute_and_extract.hpp"
#include "../sources/commute/PairAnnihilation/commute_and_extract.hpp"
#include <mrock/utility/OutputConvenience.hpp>

#include <sstream>

using namespace NickelCUT::commute;

int main(int /*argc*/, char** /*argv*/) {
    std::ostringstream oss;

    Hamiltonian::commute_and_extract(oss);
    //PairAnnihilation::commute_and_extract(oss);
    
    //////////////////////////////////////////////////////////////

    mrock::utility::save_string_raw(oss.str(), "commute_output.tex");

    return 0;
}