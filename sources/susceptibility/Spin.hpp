#pragma once

#include "ieom_base.hpp"

namespace NickelCUT::susceptibility {

class Spin : public ieom_base {
public:
    Spin(const FlowContainer& _flow_result, const Model& model);
    Spin(const ieom_base* const other_instance);

    void fill_chi_matrix(const momentum& x) final;
    void fill_U_matrix(const momentum& x) final;
};

}  // namespace NickelCUT::susceptibility
