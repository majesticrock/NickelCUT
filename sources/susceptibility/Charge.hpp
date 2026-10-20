#pragma once

#include "ieom_base.hpp"

namespace NickelCUT::susceptibility {

class Charge : public ieom_base {
public:
    Charge(const flow::FlowContainer& _flow_result, const flow::Model& model);
    Charge(const ieom_base* const other_instance);

    void fill_M(const momentum& x) final;
    void fill_N(const momentum& x) final;
};

}  // namespace NickelCUT::susceptibility
