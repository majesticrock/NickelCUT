#include "WickOrderedCollector.hpp"

#include <mrock/symbolic_operators/WickSymmetry.hpp>

#include <memory>
#include <ostream>
#include <utility>
#include <vector>
#include <list>

namespace mrock::symbolic_operators::experimental {

// Computes n!!
template <typename IntegerType, typename RealType = std::size_t>
constexpr RealType double_factorial(const IntegerType n) {
    RealType result{1};
    for (std::make_signed_t<IntegerType> i = n; i > 0; i -= 2) {
        result *= i;
    }
    return result;
}

using IndexVector = std::vector<std::size_t>;
using Pair = std::pair<std::size_t, std::size_t>;
using Pairing = std::vector<Pair>;
using PairingVector = std::vector<Pairing>;

void generate_pairings(
    const IndexVector& remaining,
    std::size_t n,
    Pairing& current,
    PairingVector& result
) {
    // We have generated enough pairs
    if (current.size() == n) {
        result.push_back(current);
        return;
    }

    // Not enough elements left to make another pair
    if (remaining.size() < 2) {
        return;
    }

    for (size_t i = 0U; i < remaining.size(); ++i) {
        std::size_t first = remaining[i];
        for (size_t j = i + 1U; j < remaining.size(); ++j) {
            std::size_t second = remaining[j];

            // Build the list of unused indices after removing first and second
            IndexVector newRemaining;
            for (size_t k = 1; k < remaining.size(); k++) {
                if (k != i && k != j) {
                    newRemaining.push_back(remaining[k]);
                }
            }

            current.push_back({first, second});

            generate_pairings(newRemaining, n, current, result);

            current.pop_back();
        }
    }
}

PairingVector get_pairings(std::size_t n_indizes, std::size_t n_pairings) {
    PairingVector result;
    Pairing current;

    IndexVector indices(n_indizes);
    std::iota(indices.begin(), indices.end(), std::size_t{});

    generate_pairings(indices, n_pairings, current, result);

    return result;
}


bool permutation_is_odd(const Pairing& pairs)
{
    std::vector<int> p;
    p.reserve(2 * pairs.size());

    for (auto [a, b] : pairs) {
        p.push_back(a);
        p.push_back(b);
    }

    bool odd = false;

    for (size_t i = 0U; i < p.size(); ++i) {
        for (size_t j = i + 1; j < p.size(); ++j) {
            if (p[i] > p[j]) {
                odd = !odd;   // toggle parity
            }
        }
    }

    return odd;
}

template <typename T>
std::vector<T> remove_elements(
    const std::vector<T>& vec,
    const Pairing& pairing
) {
    std::vector<bool> remove(vec.size(), false);

    for (auto [a, b] : pairing) {
        remove[a] = true;
        remove[b] = true;
    }

    std::vector<T> result;
    result.reserve(vec.size());

    for (size_t i = 0; i < vec.size(); i++) {
        if (!remove[i]) {
            result.push_back(vec[i]);
        }
    }

    return result;
}

void wick_decompose_internal(const std::vector<Operator>& operators, 
    const std::vector<WickOperatorTemplate>& templates,
    WickOrderedTerm&& attach_to,
    WickOrderedCollector& result)
{
    result.emplace_back(std::move(attach_to));
    result.back().wick_expression.operators = operators;

    for(std::size_t n_contractions = 1U; n_contractions <= operators.size() / 2; ++n_contractions) {
        const PairingVector pairings = get_pairings(operators.size(), n_contractions);

        for (const auto& pairing : pairings) {
            std::list<TemplateResult> template_results;

            for (const auto& pair : pairing) {
                for (const auto& wick_template : templates) {
                    template_results.emplace_back(wick_template.create_from_operators(operators[pair.first], operators[pair.second]));
                }
            }

            if(!std::all_of(template_results.begin(), template_results.end(), [](const auto& tr) { return !!tr;}))
            {
                continue;
            }

            WickOrderedTerm new_term = WickOrderedTerm::from_wick_term_skip_wick_expression(attach_to);
            new_term.multiplicity *= permutation_is_odd(pairing) ? -1 : 1;

            for (const auto& tr : template_results) {
                new_term.delta_momenta.emplace_back(tr.momentum_delta);
                new_term.include_template_result(tr.results.front());
            }

            new_term.wick_expression.operators = remove_elements(operators, pairing);
            result.emplace_back(std::move(new_term));
        }
    }
}

WickOrderedCollector wick_decompose(const Term& term, const std::vector<WickOperatorTemplate>& templates)
{
    const std::size_t estimated_size = double_factorial(term.get_operators().size()) * templates.size();
    WickOrderedCollector result;
    result.reserve(estimated_size);

    WickOrderedTerm initial(term);

    wick_decompose_internal(term.get_operators(), templates, std::move(initial), result);

    return result;
}

void clean_wick_ordered_terms(WickOrderedCollector& terms,
    const std::vector<std::unique_ptr<WickSymmetry>>& symmetries /*= std::vector<std::unique_ptr<WickSymmetry>>{}*/)
{
    for (auto& term : terms) {
        for (std::vector<Coefficient>::iterator it = term.coefficients.begin(); it != term.coefficients.end();) {
            if (it->name == "") {
                it = term.coefficients.erase(it);
            } else {
                ++it;
            }
        }
    }

    for (WickOrderedCollector::iterator it = terms.begin(); it != terms.end();) {
        if (!(it->resolve_deltas())) {
            it = terms.erase(it);
            continue;
        }
        it->discard_zero_momenta();
        it->rename_sums();
        it->sort();

        // TODO: Maybe is_pauli_forbidden()

        for (const auto& symmetry : symmetries) {
            symmetry->apply_to(it->operators);
        }

        for (auto& coeff : it->coefficients) {
            coeff.apply_custom_symmetry();
        }
        ++it;
    }

    // Setup so that we always have a structure like delta_(l,k+something)
    for (auto& term : terms) {
        for (auto& delta : term.delta_momenta) {
            assert(delta.first.momentum_list.size() == 1U);
            int l_is = delta.first.is_used_at('l');
            if (l_is == 0)
                continue;

            l_is = delta.second.is_used_at('l');
            if (l_is == -1) {
                std::cout << "################\n# Broken term\n";
                std::cout << term << std::endl;
                throw std::runtime_error("There is no l in the delta, but we need an l!");
            }
            const Momentum l_mom('l', delta.second.momentum_list[l_is].factor);
            const Momentum remainder = delta.second - l_mom;
            delta -= remainder;
            std::swap(delta.first, delta.second);
            if (delta.first.add_Q) {
                delta.second.add_Q = !delta.second.add_Q;
                delta.first.add_Q = false;
            }
        }
    }
    // remove duplicates
    for (std::size_t i = 0U; i < terms.size(); i++) {
        for (std::size_t j = i + 1U; j < terms.size(); j++) {
            if (terms[i] == terms[j]) {
                terms[i].multiplicity += terms[j].multiplicity;
                terms.erase(terms.begin() + j);
                --i;
                break;
            }
        }
    }
    // removes any terms that have a 0 prefactor
    for (auto it = terms.begin(); it != terms.end();) {
        if (it->multiplicity == 0) {
            it = terms.erase(it);
        } else {
            ++it;
        }
    }

    auto predicate = [](const WickOrderedTerm& left, const WickOrderedTerm& right) -> bool {
        if (left.delta_momenta.empty() && right.delta_momenta.size() > 0) {
            return true;
        } else if (left.delta_momenta.size() > 0 && right.delta_momenta.size() > 0) {
            if (left.delta_momenta.size() < right.delta_momenta.size()) {
                return true;
            } else if (left.delta_momenta.size() == right.delta_momenta.size()) {
                if (left.delta_momenta[0].second.add_Q && !(right.delta_momenta[0].second.add_Q)) {
                    return true;
                } else if (!left.coefficients.empty() && right.coefficients[0].name < left.coefficients[0].name) {
                    return true;
                } else if ((!left.coefficients.empty() && right.coefficients[0].name == left.coefficients[0].name) ||
                           left.coefficients.empty()) {
                    if (!left.operators.empty() && right.operators.empty()) {
                        return true;
                    } else if ((!left.operators.empty() && !right.operators.empty()) &&
                               left.operators.front().type < right.operators.front().type) {
                        return true;
                    }
                }
            }
        } else if (left.delta_momenta.empty() && right.delta_momenta.empty()) {
            if (left.coefficients.size() < right.coefficients.size()) {
                return true;
            } else if (!left.coefficients.empty() && !right.coefficients.empty()) {
                if (right.coefficients[0].name < left.coefficients[0].name) {
                    return true;
                } else if (right.coefficients[0].name == left.coefficients[0].name) {
                    if (left.operators.size() > right.operators.size()) {
                        return true;
                    } else if ((!left.operators.empty() && !right.operators.empty()) &&
                               left.operators.front().type < right.operators.front().type) {
                        return true;
                    }
                }
            }
        }
        return false;
    };

    // Sort terms
    for (std::size_t i = 0; i < terms.size(); i++) {
        for (std::size_t j = i + 1; j < terms.size(); j++) {
            if (predicate(terms[i], terms[j]))
                std::swap(terms[i], terms[j]);
        }
    }
}

std::ostream& operator<<(std::ostream& os, const WickOrderedCollector& terms) {
    for (WickOrderedCollector::const_iterator it = terms.begin(); it != terms.end(); ++it) {
        os << "\t&" << *it;
        if (it != terms.end() - 1) {
            os << " \\\\";
        }
        os << "\n";
    }
    return os;
}
}