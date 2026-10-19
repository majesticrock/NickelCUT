#pragma once

#include "ieom_base.hpp"

namespace NickelCUT::susceptibility
{

class Spin : public ieom_base {
public:
    Spin(const flow::FlowContainer& _flow_result, const flow::Model& model);
    Spin(const ieom_base * const other_instance);

    void fill_M(const momentum& x) final;
    void fill_N(const momentum& x) final;
};

} // namespace NickelCUT::susceptibility
