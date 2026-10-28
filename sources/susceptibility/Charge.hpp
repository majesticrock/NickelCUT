#pragma once

#include "ieom_base.hpp"

namespace NickelCUT::susceptibility {

class Charge : public ieom_base {
public:
    Charge(const FlowContainer& _flow_result, const Model& model);
    Charge(const ieom_base* const other_instance);

    void fill_chi_matrix(const momentum& x) final;
    void fill_U_matrix(const momentum& x) final;
};

}  // namespace NickelCUT::susceptibility
