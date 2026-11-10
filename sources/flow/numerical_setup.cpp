#include "numerical_setup.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace NickelCUT::flow {
double finite_size_energy_gap(const Model& model) {
    constexpr int fine_mesh_size = 201;
    constexpr double energy_tolerance = 1e-12;
    constexpr double pi = std::numbers::pi;
    constexpr double two_pi = 2. * pi;
    constexpr double maximum_contour_distance_factor = 0.99 * pi * std::numbers::sqrt2;

    std::array<double, fine_mesh_size> fine_momenta{};
    for (int i = 0; i < fine_mesh_size; ++i) {
        fine_momenta[i] = -pi + two_pi * i / (fine_mesh_size - 1);
    }

    const auto energy_at = [&model](double kx, double ky) {
        return -2. * (std::cos(kx) + std::cos(ky) + 2. * model.tprime * std::cos(kx) * std::cos(ky));
    };
    std::vector<double> fine_energies(fine_mesh_size * fine_mesh_size);
    double minimum_energy = std::numeric_limits<double>::infinity();
    double maximum_energy = -std::numeric_limits<double>::infinity();
    for (int y = 0; y < fine_mesh_size; ++y) {
        for (int x = 0; x < fine_mesh_size; ++x) {
            const double energy = energy_at(fine_momenta[x], fine_momenta[y]);
            fine_energies[y * fine_mesh_size + x] = energy;
            minimum_energy = std::min(minimum_energy, energy);
            maximum_energy = std::max(maximum_energy, energy);
        }
    }

    if (model.E_F < minimum_energy || model.E_F > maximum_energy) {
        throw std::invalid_argument("Fermi level is outside the dispersion bandwidth.");
    }

    const auto fine_energy = [&fine_energies](int x, int y) {
        return fine_energies[y * fine_mesh_size + x];
    };
    using Momentum = std::pair<double, double>;
    std::vector<Momentum> fermi_surface_points;
    const auto add_edge_intersection = [&fermi_surface_points, &fine_momenta, &fine_energy, &model](
                                          int x1, int y1, int x2, int y2) {
        const double energy1 = fine_energy(x1, y1);
        const double energy2 = fine_energy(x2, y2);
        if ((energy1 < model.E_F && energy2 >= model.E_F) ||
            (energy2 < model.E_F && energy1 >= model.E_F)) {
            const double fraction = (model.E_F - energy1) / (energy2 - energy1);
            fermi_surface_points.emplace_back(
                fine_momenta[x1] + fraction * (fine_momenta[x2] - fine_momenta[x1]),
                fine_momenta[y1] + fraction * (fine_momenta[y2] - fine_momenta[y1]));
        }
    };

    for (int y = 0; y < fine_mesh_size - 1; ++y) {
        for (int x = 0; x < fine_mesh_size - 1; ++x) {
            add_edge_intersection(x, y, x + 1, y);
            add_edge_intersection(x + 1, y, x + 1, y + 1);
            add_edge_intersection(x + 1, y + 1, x, y + 1);
            add_edge_intersection(x, y + 1, x, y);
        }
    }

    if (fermi_surface_points.empty()) {
        double closest_energy_distance = std::numeric_limits<double>::infinity();
        int closest_x = 0;
        int closest_y = 0;
        for (int y = 0; y < fine_mesh_size; ++y) {
            for (int x = 0; x < fine_mesh_size; ++x) {
                const double distance = std::abs(fine_energy(x, y) - model.E_F);
                if (distance <= energy_tolerance) {
                    fermi_surface_points.emplace_back(fine_momenta[x], fine_momenta[y]);
                }
                if (distance < closest_energy_distance) {
                    closest_energy_distance = distance;
                    closest_x = x;
                    closest_y = y;
                }
            }
        }
        if (fermi_surface_points.empty()) {
            fermi_surface_points.emplace_back(fine_momenta[closest_x], fine_momenta[closest_y]);
        }
    }

    const double maximum_contour_distance = maximum_contour_distance_factor / L;
    double energy_gap = 0.;
    bool has_fermi_surface_momentum = false;
    for (int y = 0; y < L; ++y) {
        for (int x = 0; x < L; ++x) {
            const double kx = -pi + two_pi * x / L;
            const double ky = -pi + two_pi * y / L;
            double minimum_distance_squared = std::numeric_limits<double>::infinity();
            for (const auto& [contour_kx, contour_ky] : fermi_surface_points) {
                const double delta_kx = std::remainder(kx - contour_kx, two_pi);
                const double delta_ky = std::remainder(ky - contour_ky, two_pi);
                minimum_distance_squared = std::min(
                    minimum_distance_squared, delta_kx * delta_kx + delta_ky * delta_ky);
            }
            if (std::sqrt(minimum_distance_squared) > maximum_contour_distance) {
                continue;
            }

            has_fermi_surface_momentum = true;
            const double energy = energy_at(kx, ky);
            const std::array<std::pair<int, int>, 4> neighbors{{
                {(x + L - 1) % L, y},
                {(x + 1) % L, y},
                {x, (y + L - 1) % L},
                {x, (y + 1) % L},
            }};
            for (const auto& [neighbor_x, neighbor_y] : neighbors) {
                const double neighbor_kx = -pi + two_pi * neighbor_x / L;
                const double neighbor_ky = -pi + two_pi * neighbor_y / L;
                energy_gap = std::max(energy_gap, std::abs(energy - energy_at(neighbor_kx, neighbor_ky)));
            }
        }
    }

    if (!has_fermi_surface_momentum || !std::isfinite(energy_gap) || energy_gap <= 0.) {
        throw std::runtime_error("Could not determine a positive finite-size energy gap at the Fermi surface.");
    }
    return energy_gap;
}
}  // namespace NickelCUT::flow
