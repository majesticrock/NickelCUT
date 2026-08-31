#include "extract_flow_coefficients.hpp"

#include "../../experimental/WickOrderedCollector.hpp"

#include <mrock/symbolic_operators/WickTermCollector.hpp>

#include <array>

namespace NickelCUT::commute::Hamiltonian {
using namespace mrock::symbolic_operators;

std::array<WickTermCollector, 3> extract_flow_coefficients(const experimental::WickOrderedCollector& terms) {
    std::array<WickTermCollector, 3> flow_coefficients;
    flow_coefficients[0].reserve(std::count_if(terms.begin(), terms.end(), 
        [](const experimental::WickOrderedTerm& term){
            return term.is_bilinear();
    }));
    flow_coefficients[1].reserve(std::count_if(terms.begin(), terms.end(), 
        [](const experimental::WickOrderedTerm& term){
            return term.is_bilinear();
    }));
    flow_coefficients[2].reserve(std::count_if(terms.begin(), terms.end(), 
        [](const experimental::WickOrderedTerm& term){
            return term.is_quartic();
    }));

    auto wick_term_from_normal_ordered = [](const experimental::WickOrderedTerm& term) {
        return WickTerm(
                term.multiplicity,
                term.coefficients,
                term.sums,
                term.delta_momenta,
                term.delta_indices,
                term.operators
            );
    };

    const std::vector<Index> removeable_indices_bilinear = {Index::Sigma};
    const std::vector<Index> removeable_indices_quartic = {Index::Sigma, Index::SigmaPrime};

    const std::vector<MomentumSymbol::name_type> removeable_momenta_bilinear = {'K'};
    const std::vector<MomentumSymbol::name_type> removeable_momenta_quartic = {'K', 'P', 'Q'};

    auto remove_sums = [](WickTerm& term, const std::vector<Index>& indices, const std::vector<MomentumSymbol::name_type>& momenta) {
        std::erase_if(term.sums.spins.summations, [&indices](Index index) {
            return exists_in(indices, index);
        });

        std::erase_if(term.sums.momenta.summations, [&momenta](MomentumSymbol::name_type momentum) {
            return exists_in(momenta, momentum);
        });
    };

    for (const auto& term : terms) {
        if(term.is_bilinear()) {
            flow_coefficients[0].push_back(wick_term_from_normal_ordered(term));
            remove_sums(flow_coefficients[0].back(), removeable_indices_bilinear, removeable_momenta_bilinear);
        }
        else if (term.is_quartic()) {
            flow_coefficients[1].push_back(wick_term_from_normal_ordered(term));
            remove_sums(flow_coefficients[1].back(), removeable_indices_quartic, removeable_momenta_quartic);
        }
    }

    flow_coefficients[2] = flow_coefficients[1];
    for (auto& term : flow_coefficients[1]) {
        term.replace_each_index(Index::Sigma, Index::SpinUp);
        term.replace_each_index(Index::SigmaPrime, Index::SpinDown);
    }

    for (auto& term : flow_coefficients[2]) {
        term.replace_each_index(Index::Sigma, Index::SpinUp);
        term.replace_each_index(Index::SigmaPrime, Index::SpinUp);
    }

    return flow_coefficients;
} 

} // namespace NickelCUT::commute