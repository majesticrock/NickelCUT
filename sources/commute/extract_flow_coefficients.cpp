#include "extract_flow_coefficients.hpp"

#include "../experimental/WickOrderedCollector.hpp"

#include <mrock/symbolic_operators/WickTerm.hpp>

#include <array>

namespace NickelCUT::commute {
using namespace mrock::symbolic_operators;

std::array<WickTermCollector, 3> extract_flow_coefficients(const experimental::WickOrderedCollector& terms) {
    std::array<WickTermCollector, 3> flow_coefficients;
    flow_coefficients[0].reserve(std::count_if(terms.begin(), terms.end(), 
        [](const experimental::WickOrderedTerm& term){
            return term.wick_expression.empty();
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
                term.delta_indizes,
                term.operators
            );
    };

    const std::vector<Index> removeable_indizes_bilinear = {Index::Sigma};
    const std::vector<Index> removeable_indizes_quartic = {Index::Sigma, Index::SigmaPrime};

    const std::vector<MomentumSymbol::name_type> removeable_momenta_bilinear = {'p'};
    const std::vector<MomentumSymbol::name_type> removeable_momenta_quartic = {'p', 'q', 'r'};

    auto remove_sums = [](WickTerm& term, const std::vector<Index>& indizes, const std::vector<MomentumSymbol::name_type>& momenta) {
        std::erase_if(term.sums.spins.summations, [&indizes](Index index) {
            return exists_in(indizes, index);
        });

        std::erase_if(term.sums.momenta.summations, [&momenta](MomentumSymbol::name_type momentum) {
            return exists_in(momenta, momentum);
        });
    };

    for (const auto& term : terms) {
        if (term.wick_expression.empty()) {
            flow_coefficients[0].push_back(wick_term_from_normal_ordered(term));
        }
        else if(term.is_bilinear()) {
            flow_coefficients[1].push_back(wick_term_from_normal_ordered(term));
            remove_sums(flow_coefficients[1].back(), removeable_indizes_bilinear, removeable_momenta_bilinear);
        }
        else if (term.is_quartic()) {
            flow_coefficients[2].push_back(wick_term_from_normal_ordered(term));
            remove_sums(flow_coefficients[2].back(), removeable_indizes_quartic, removeable_momenta_quartic);
        }
    }

    return flow_coefficients;
} 

} // namespace NickelCUT::commute