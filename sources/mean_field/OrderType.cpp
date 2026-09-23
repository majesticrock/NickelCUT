#include "OrderType.hpp"

namespace NickelCUT::mean_field
{

std::string to_string(OrderType order)
{
    const auto value = static_cast<unsigned char>(order);

    if (value == static_cast<unsigned char>(OrderType::NormalMetal)) {
        return "NM";
    }
    if (value == static_cast<unsigned char>(OrderType::UnrecognizedOrder)) {
        return "UNRECOGNIZED";
    }

    std::string result;
    auto append = [&result](const char* name) {
        if (!result.empty()) {
            result += "+";
        }
        result += name;
    };

    if ((value & static_cast<unsigned char>(OrderType::swave_SC)) != 0) {
        append("s-SC");
    }
    if ((value & static_cast<unsigned char>(OrderType::dwave_SC)) != 0) {
        append("d-SC");
    }
    if ((value & static_cast<unsigned char>(OrderType::AFM)) != 0) {
        append("AFM");
    }
    if ((value & static_cast<unsigned char>(OrderType::CDW)) != 0) {
        append("CDW");
    }

    return result.empty() ? "UNRECOGNIZED" : result;
}

} // namespace NickelCUT::mean_field
