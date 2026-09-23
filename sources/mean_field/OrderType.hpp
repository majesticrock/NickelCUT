#pragma once

#include "../flow/Model.hpp"
#include "../flow/ExtractionContainer.hpp"
#include "../flow/momentum_iterator.hpp"
#include "ModelAttributes.hpp"
#include "../L.hpp"

#include <Eigen/Dense>

#include <cassert>
#include <string>

namespace NickelCUT::mean_field {

enum class OrderType : unsigned char {
    NormalMetal = 0,
    swave_SC = 1,
    dwave_SC = 2,
    AFM = 4,
    CDW = 8,
    UnrecognizedOrder = 255
};

std::string to_string(OrderType order);

constexpr OrderType operator|(OrderType lhs, OrderType rhs) 
{
    return static_cast<OrderType>(static_cast<unsigned char>(lhs) | static_cast<unsigned char>(rhs));
}
constexpr OrderType operator&(OrderType lhs, OrderType rhs) 
{
    return static_cast<OrderType>(static_cast<unsigned char>(lhs) & static_cast<unsigned char>(rhs));
}

} // namespace NickelCUT::mean_field