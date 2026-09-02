#include "BookKeeper.hpp"
#include "DecouplingChannel.hpp"
#include "FlowContainer.hpp"

#include <mrock/utility/OutputConvenience.hpp>

#include <utility>
#include <deque>
#include <iostream>

namespace NickelCUT::flow
{

ExtractionContainer::ExtractionContainer(const FlowContainer& x)
    : single_particle_energy(DecouplingChannel::SingleParticleEnergy(x)),
    density_wave(DecouplingChannel::DensityWave(x)),
    superconductivity(DecouplingChannel::Superconductivity(x)),
    dispersion(x.dispersion)
{ }

BookKeeper::BookKeeper(const FlowContainer& initial_flow_state, double _dl) 
    : lowest_ROD{ initial_flow_state.residual_offdiagonality() },
    l_of_lowest_ROD{ 0.0 },
    index_of_lowest_ROD{ 0U },
    l_times{ 0.0 },
    residual_offdiagonalities{ lowest_ROD },
    extracted_channels{ ExtractionContainer(initial_flow_state) },
    lowest_ROD_state{ initial_flow_state },
    dl{ _dl },
    min_ROD_difference{ 0.02 * lowest_ROD },
    begin(clock::now()), 
    last(begin),
    current_idx{ 0U }
{
    std::cout << mrock::utility::time_stamp() << "   -   " << "Starting calculations...\n"
        << "Initial ROD = " << lowest_ROD << "        spacing saves at at least dl = " << dl << std::endl;
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

void BookKeeper::operator()(const FlowContainer &x, double l)
{
    if (l - l_times.back() < dl) return;

    const double current_ROD = x.residual_offdiagonality();

    if (process_step(l, residual_offdiagonalities.back())) {
        lowest_ROD_state = x;

        if (std::abs(current_ROD - residual_offdiagonalities.back()) < min_ROD_difference) {
            l_times.back() = l;
            residual_offdiagonalities.back() = current_ROD;
            extracted_channels.back() = ExtractionContainer(x);
        }
    }
    else if (std::abs(current_ROD - residual_offdiagonalities.back()) > min_ROD_difference) {
        l_times.push_back(l);
        residual_offdiagonalities.push_back(current_ROD);
        extracted_channels.push_back(ExtractionContainer(x));
    }

    if (residual_offdiagonalities.back() > 5 * residual_offdiagonalities.front()) {
        throw LargeRODException();
    }
}

void to_json(nlohmann::json& j, const ExtractionContainer& extracted_channels) noexcept
{
    j = nlohmann::json{
        { "single_particle_energy_differing",   extracted_channels.single_particle_energy.first.as_2D_array()  },
        { "single_particle_energy_same",        extracted_channels.single_particle_energy.second.as_2D_array() },
        { "density_wave_differing",             extracted_channels.density_wave.first.as_2D_array()            },
        { "density_wave_same",                  extracted_channels.density_wave.second.as_2D_array()           },
        { "superconductivity",                  extracted_channels.superconductivity.as_2D_array()             },
        { "dispersion",                         extracted_channels.dispersion                                  }
    };
}

void to_json(nlohmann::json& j, const BookKeeper& book_keeper) noexcept
{
    j = nlohmann::json{
        { "l_times",                    book_keeper.l_times                   },
        { "number_of_data_points",      book_keeper.l_times.size()            },
        { "residual_offdiagonalities",  book_keeper.residual_offdiagonalities },
        { "lowest_ROD",                 book_keeper.lowest_ROD                },
        { "index_of_lowest_ROD",        book_keeper.index_of_lowest_ROD       },
        { "l_of_lowest_ROD",            book_keeper.l_of_lowest_ROD           },
        { "extracted_channels",         book_keeper.extracted_channels        },
        { "lowest_ROD_state",           book_keeper.lowest_ROD_state          }
    };
}

} // namespace NickelCUT::flow
