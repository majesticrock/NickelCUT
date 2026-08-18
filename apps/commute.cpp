#include "../sources/commute/Hamiltonian/commute_and_print.hpp"
#include "../sources/commute/PairAnnihilation/commute_and_print.hpp"
#include <mrock/utility/OutputConvenience.hpp>

#include <sstream>

using namespace NickelCUT::commute;

int main(int /*argc*/, char** /*argv*/) {
    std::ostringstream oss;

    Hamiltonian::commute_and_print(oss);
    PairAnnihilation::commute_and_print(oss);
    
    //////////////////////////////////////////////////////////////

    mrock::utility::save_string_raw(oss.str(), "commute_output.tex");

    return 0;
}