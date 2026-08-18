#include "BookKeeper.hpp"

#include <mrock/utility/OutputConvenience.hpp>

#include <iostream>

namespace NickelCUT::flow
{

BookKeeper::BookKeeper(double initial_ROD) 
    : lowest_ROD{ initial_ROD },
    l_of_lowest_ROD{ 0.0 },
    index_of_lowest_ROD{ 0U },
    begin(clock::now()), 
    last(begin),
    current_idx{ 0U }
{
    std::cout << mrock::utility::time_stamp() << "   -   " << "Starting calculations...\n"
        << "Initial ROD = " << initial_ROD << std::endl;
};

bool BookKeeper::process_step(double current_l, double ROD) {
    bool updated = false;
    ++current_idx;
    if (ROD < lowest_ROD) {
        lowest_ROD = ROD;
        l_of_lowest_ROD = current_l;
        index_of_lowest_ROD = current_idx;
        updated = true;
    }
    clock::time_point now = clock::now();
    std::cout << "//------------------------------------------------------//\n"
        << "Step #" << current_idx << "\t" << mrock::utility::time_stamp() << "\n"
        << "l = " << current_l
        << "\t\tROD = " << ROD << "\n"
        << "Step took " << std::chrono::duration_cast<std::chrono::milliseconds>(now - last).count() << "ms to execute."
        << std::endl;
    last = now;
    return updated;
}

void BookKeeper::print_final() const {
    clock::time_point now = clock::now();
    std::cout << "//------------------------------------------------------//\n"
        << "\t Flow program finished at "
        << mrock::utility::time_stamp() << "\n"
        << "lowest ROD achieved after " << index_of_lowest_ROD << " steps at l=" << l_of_lowest_ROD << "."
        << "\t\tlowest ROD = " << lowest_ROD << "\n"
        << "Total executation took " << std::chrono::duration_cast<std::chrono::seconds>(now - begin).count() << "s.\n"
        << "Goodbye."
        << std::endl;
}

} // namespace NickelCUT::flow
