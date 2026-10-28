#pragma once

#include "ieom_base.hpp"

namespace NickelCUT::susceptibility {

class SC : public ieom_base {
public:
    SC(const FlowContainer& _flow_result, const Model& model);
    SC(const ieom_base* const other_instance);

    void fill_chi_matrix(const momentum& x) final;
    void fill_U_matrix(const momentum& x) final;
};

}  // namespace NickelCUT::susceptibility
