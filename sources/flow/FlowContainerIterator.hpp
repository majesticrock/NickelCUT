#pragma once

#include <cstddef>
#include <iterator>
#include <type_traits>

namespace NickelCUT::flow {

template<typename ContainerT>
class FlowContainerIterator {
public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = double;
    using difference_type = std::ptrdiff_t;
    using reference = std::conditional_t<std::is_const_v<ContainerT>, const double&, double&>;
    using pointer = std::conditional_t<std::is_const_v<ContainerT>, const double*, double*>;

    constexpr FlowContainerIterator() noexcept = default;
    constexpr FlowContainerIterator(ContainerT* container, std::size_t index) noexcept
        : _container(container), _index(index) {}

    constexpr reference operator*() const noexcept {
        return _container->value_at(_index);
    }

    constexpr pointer operator->() const noexcept {
        return std::addressof(_container->value_at(_index));
    }

    constexpr FlowContainerIterator& operator++() noexcept {
        ++_index;
        return *this;
    }

    constexpr FlowContainerIterator operator++(int) noexcept {
        auto copy = *this;
        ++*this;
        return copy;
    }

    constexpr FlowContainerIterator& operator--() noexcept {
        --_index;
        return *this;
    }

    constexpr FlowContainerIterator operator--(int) noexcept {
        auto copy = *this;
        --*this;
        return copy;
    }

    constexpr FlowContainerIterator& operator+=(difference_type offset) noexcept {
        _index += static_cast<std::size_t>(offset);
        return *this;
    }

    constexpr FlowContainerIterator operator+(difference_type offset) const noexcept {
        return FlowContainerIterator(_container, _index + static_cast<std::size_t>(offset));
    }

    constexpr FlowContainerIterator& operator-=(difference_type offset) noexcept {
        _index -= static_cast<std::size_t>(offset);
        return *this;
    }

    constexpr FlowContainerIterator operator-(difference_type offset) const noexcept {
        return FlowContainerIterator(_container, _index - static_cast<std::size_t>(offset));
    }

    constexpr reference operator[](difference_type offset) const noexcept {
        return _container->value_at(_index + static_cast<std::size_t>(offset));
    }

    template<typename OtherContainerT>
    friend constexpr bool operator==(const FlowContainerIterator& lhs, const FlowContainerIterator<OtherContainerT>& rhs) noexcept {
        return static_cast<const void*>(lhs._container) == static_cast<const void*>(rhs._container) && lhs._index == rhs._index;
    }

    template<typename OtherContainerT>
    friend constexpr bool operator!=(const FlowContainerIterator& lhs, const FlowContainerIterator<OtherContainerT>& rhs) noexcept {
        return !(lhs == rhs);
    }

private:
    ContainerT* _container = nullptr;
    std::size_t _index = 0;
};

template<typename ContainerT>
constexpr FlowContainerIterator<ContainerT> operator+(typename FlowContainerIterator<ContainerT>::difference_type offset,
                                                       FlowContainerIterator<ContainerT> iterator) noexcept {
    return iterator += offset;
}

} // namespace NickelCUT::flow
