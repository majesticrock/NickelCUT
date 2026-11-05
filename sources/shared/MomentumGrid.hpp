#pragma once

#include "../L.hpp"
#include "../momentum_iterator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace NickelCUT {

struct WeightedMomentum {
    momentum_iterator<L> momentum;
    double weight;
};

class MomentumGrid {
    struct State {
        bool initialized{};
        std::uint64_t generation{};
        std::vector<WeightedMomentum> fermions;
        std::vector<WeightedMomentum> transfers;
        std::array<int, N> nearest_fermion{};
        std::array<int, N> nearest_transfer{};
        std::array<bool, N> is_fermion{};
        std::array<bool, N> is_transfer{};
    };

    static State& state() {
        static State value;
        return value;
    }

    static double dispersion(double kx, double ky, double tprime) noexcept {
        return -2.0 * (std::cos(kx) + std::cos(ky) + 2.0 * tprime * std::cos(kx) * std::cos(ky));
    }

    static double periodic_distance_squared(double x1, double y1, double x2, double y2) noexcept {
        const double dx = std::remainder(x1 - x2, 2.0 * std::numbers::pi);
        const double dy = std::remainder(y1 - y2, 2.0 * std::numbers::pi);
        return dx * dx + dy * dy;
    }

    static std::vector<std::pair<double, double>> fermi_contour(double tprime, double chemical_potential) {
        constexpr int mesh_size = 201;
        constexpr int cells = mesh_size - 1;
        const double step = 2.0 * std::numbers::pi / cells;
        std::array<double, mesh_size> coordinates{};
        for (int i = 0; i < mesh_size; ++i) {
            coordinates[i] = -std::numbers::pi + step * i;
        }

        std::vector<std::pair<double, double>> points;
        points.reserve(4 * cells);
        const auto energy = [=](int x, int y) {
            return dispersion(coordinates[x], coordinates[y], tprime) - chemical_potential;
        };
        const auto add_crossing = [&points](double x1, double y1, double e1, double x2, double y2, double e2) {
            if (e1 == 0.0) {
                points.emplace_back(x1, y1);
            } else if (e2 == 0.0) {
                points.emplace_back(x2, y2);
            } else if ((e1 < 0.0) != (e2 < 0.0)) {
                const double fraction = e1 / (e1 - e2);
                points.emplace_back(x1 + fraction * (x2 - x1), y1 + fraction * (y2 - y1));
            }
        };

        for (int y = 0; y < cells; ++y) {
            for (int x = 0; x < cells; ++x) {
                const double e00 = energy(x, y);
                const double e10 = energy(x + 1, y);
                const double e01 = energy(x, y + 1);
                const double e11 = energy(x + 1, y + 1);
                add_crossing(coordinates[x], coordinates[y], e00, coordinates[x + 1], coordinates[y], e10);
                add_crossing(coordinates[x + 1], coordinates[y], e10, coordinates[x + 1], coordinates[y + 1], e11);
                add_crossing(coordinates[x + 1], coordinates[y + 1], e11, coordinates[x], coordinates[y + 1], e01);
                add_crossing(coordinates[x], coordinates[y + 1], e01, coordinates[x], coordinates[y], e00);
            }
        }

        if (points.empty()) {
            double minimum_distance = std::numeric_limits<double>::infinity();
            std::pair<double, double> closest{};
            for (int y = 0; y < mesh_size; ++y) {
                for (int x = 0; x < mesh_size; ++x) {
                    const double distance = std::abs(energy(x, y));
                    if (distance < minimum_distance) {
                        minimum_distance = distance;
                        closest = {coordinates[x], coordinates[y]};
                    }
                }
            }
            points.push_back(closest);
        }
        return points;
    }

    static bool near_contour(int x,
                             int y,
                             const std::vector<std::pair<double, double>>& contour,
                             double maximum_distance) noexcept {
        const auto& momentum = momentum_iterator<L>::momentum_cache;
        const double kx = momentum[x];
        const double ky = momentum[y];
        const double maximum_distance_squared = maximum_distance * maximum_distance;
        for (const auto& point : contour) {
            if (periodic_distance_squared(kx, ky, point.first, point.second) <= maximum_distance_squared) {
                return true;
            }
        }
        return false;
    }

    static void fill_nearest_lookup(const std::vector<WeightedMomentum>& points,
                                    std::array<int, N>& lookup) {
        for (int position = 0; position < N; ++position) {
            const momentum_iterator<L> target(position);
            int best_distance = std::numeric_limits<int>::max();
            int nearest = 0;
            for (const auto& point : points) {
                const int raw_dx = std::abs(target._x - point.momentum._x);
                const int raw_dy = std::abs(target._y - point.momentum._y);
                const int dx = std::min(raw_dx, L - raw_dx);
                const int dy = std::min(raw_dy, L - raw_dy);
                const int distance = dx * dx + dy * dy;
                if (distance < best_distance) {
                    best_distance = distance;
                    nearest = point.momentum.get_position();
                }
            }
            lookup[position] = nearest;
        }
    }

public:
    static void reset() noexcept {
        State next;
        next.generation = state().generation + 1;
        state() = std::move(next);
    }

    static void initialize(double tprime, double chemical_potential) {
        if (L < 2 || L % 2 != 0) {
            throw std::invalid_argument("The adaptive momentum grid requires an even L >= 2.");
        }

        State next;
        next.initialized = true;
        next.generation = state().generation + 1;
        const int coarse_size = L / 2;
        const double maximum_distance = 0.99 * std::numbers::pi * std::sqrt(2.0) / coarse_size;
        const auto contour = fermi_contour(tprime, chemical_potential);
        std::array<bool, N> near_surface{};

        for (int position = 0; position < N; ++position) {
            const momentum_iterator<L> momentum(position);
            near_surface[position] = near_contour(momentum._x, momentum._y, contour, maximum_distance);
        }

        for (int position = 0; position < N; ++position) {
            const momentum_iterator<L> momentum(position);
            const int reflected_x = (L - momentum._x) % L;
            const int reflected_y = (L - momentum._y) % L;
            // Keep the selected set closed under the symmetries used by the data frames.
            const bool near_fermi_surface =
                near_surface[position] || near_surface[reflected_x + L * momentum._y] ||
                near_surface[momentum._x + L * reflected_y] || near_surface[reflected_x + L * reflected_y] ||
                near_surface[momentum._y + L * momentum._x] ||
                near_surface[reflected_y + L * momentum._x] ||
                near_surface[momentum._y + L * reflected_x] ||
                near_surface[reflected_y + L * reflected_x];
            if (near_fermi_surface) {
                next.fermions.push_back({momentum, 1.0});
                next.is_fermion[position] = true;
            } else if (momentum._x % 2 == 0 && momentum._y % 2 == 0) {
                next.fermions.push_back({momentum, 4.0});
                next.is_fermion[position] = true;
            }

            if (momentum._x % 2 == 0 && momentum._y % 2 == 0) {
                next.transfers.push_back({momentum, 4.0});
                next.is_transfer[position] = true;
            }
        }

        if (next.fermions.empty() || next.transfers.empty()) {
            throw std::runtime_error("Failed to construct the adaptive momentum grid.");
        }
        double total_fermion_weight{};
        for (const auto& point : next.fermions) {
            total_fermion_weight += point.weight;
        }
        // A coarse node represents four fine-grid cells; normalize to retain the original sum scale.
        const double weight_normalization = static_cast<double>(N) / total_fermion_weight;
        for (auto& point : next.fermions) {
            point.weight *= weight_normalization;
        }
        fill_nearest_lookup(next.fermions, next.nearest_fermion);
        fill_nearest_lookup(next.transfers, next.nearest_transfer);
        state() = std::move(next);
    }

    static bool initialized() noexcept { return state().initialized; }
    static std::uint64_t generation() noexcept { return state().generation; }

    static const std::vector<WeightedMomentum>& fermion_points() {
        auto& current = state();
        if (!current.initialized) {
            static const std::vector<WeightedMomentum> full_grid = [] {
                std::vector<WeightedMomentum> result;
                result.reserve(N);
                for (int position = 0; position < N; ++position) {
                    result.push_back({momentum_iterator<L>(position), 1.0});
                }
                return result;
            }();
            return full_grid;
        }
        return current.fermions;
    }

    static const std::vector<WeightedMomentum>& transfer_points() {
        auto& current = state();
        if (!current.initialized) {
            return fermion_points();
        }
        return current.transfers;
    }

    static momentum_iterator<L> nearest_fermion(momentum_iterator<L> momentum) noexcept {
        const auto& current = state();
        return current.initialized ? momentum_iterator<L>(current.nearest_fermion[momentum.get_position()]) : momentum;
    }

    static momentum_iterator<L> nearest_transfer(momentum_iterator<L> momentum) noexcept {
        const auto& current = state();
        return current.initialized ? momentum_iterator<L>(current.nearest_transfer[momentum.get_position()]) : momentum;
    }

    static bool is_fermion(int position) noexcept {
        const auto& current = state();
        return !current.initialized || current.is_fermion[position];
    }

    static bool is_transfer(int position) noexcept {
        const auto& current = state();
        return !current.initialized || current.is_transfer[position];
    }

};

}  // namespace NickelCUT
