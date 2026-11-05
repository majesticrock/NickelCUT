#include "../sources/shared/InteractionDataFrame.hpp"
#include "../sources/shared/MomentumGrid.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
using namespace NickelCUT;

namespace {
int reflect_x(int position) {
    const momentum_iterator<L> momentum(position);
    return (L - momentum._x) % L + L * momentum._y;
}

int reflect_y(int position) {
    const momentum_iterator<L> momentum(position);
    return momentum._x + L * ((L - momentum._y) % L);
}

int rotate_90(int position) {
    const momentum_iterator<L> momentum(position);
    return (L - momentum._y) % L + L * momentum._x;
}

int distance_squared(int first, int second) {
    const momentum_iterator<L> a(first);
    const momentum_iterator<L> b(second);
    const int raw_dx = std::abs(a._x - b._x);
    const int raw_dy = std::abs(a._y - b._y);
    const int dx = std::min(raw_dx, L - raw_dx);
    const int dy = std::min(raw_dy, L - raw_dy);
    return dx * dx + dy * dy;
}
}  // namespace

int main() {
    MomentumGrid::initialize(-0.3, -0.8);

    double fermion_weight{};
    double minimum_point_weight = std::numeric_limits<double>::infinity();
    double maximum_point_weight{};
    for (const auto& point : MomentumGrid::fermion_points()) {
        fermion_weight += point.weight;
        minimum_point_weight = std::min(minimum_point_weight, point.weight);
        maximum_point_weight = std::max(maximum_point_weight, point.weight);
        if (!MomentumGrid::is_fermion(point.momentum.get_position())) {
            std::cerr << "An adaptive fermion point is missing from the lookup mask\n";
            return 1;
        }
    }
    if (std::abs(maximum_point_weight / minimum_point_weight - 4.0) > 1e-12) {
        std::cerr << "Fine and coarse momenta do not have a 1:4 quadrature-weight ratio\n";
        return 1;
    }
    if (std::abs(fermion_weight - N) > 1e-10) {
        std::cerr << "Adaptive fermion weights do not preserve the full-grid normalization\n";
        return 1;
    }

    double transfer_weight{};
    if (MomentumGrid::transfer_points().size() != N / 4) {
        std::cerr << "The transfer grid is not L/2 by L/2\n";
        return 1;
    }
    for (const auto& point : MomentumGrid::transfer_points()) {
        transfer_weight += point.weight;
        const auto momentum = point.momentum;
        if (momentum._x % 2 != 0 || momentum._y % 2 != 0 || point.weight != 4.0) {
            std::cerr << "Unexpected transfer-grid point or weight\n";
            return 1;
        }
    }
    if (std::abs(transfer_weight - N) > 1e-10) {
        std::cerr << "Transfer-grid weights do not preserve the full-grid normalization\n";
        return 1;
    }

    for (int position = 0; position < N; ++position) {
        if (MomentumGrid::is_fermion(position) != MomentumGrid::is_fermion(reflect_x(position)) ||
            MomentumGrid::is_fermion(position) != MomentumGrid::is_fermion(reflect_y(position)) ||
            MomentumGrid::is_fermion(position) != MomentumGrid::is_fermion(rotate_90(position))) {
            std::cerr << "The adaptive fermion grid is not square-lattice symmetric\n";
            return 1;
        }

        const auto nearest_fermion = MomentumGrid::nearest_fermion(momentum_iterator<L>(position));
        const auto nearest_transfer = MomentumGrid::nearest_transfer(momentum_iterator<L>(position));
        if (!MomentumGrid::is_fermion(nearest_fermion.get_position()) ||
            !MomentumGrid::is_transfer(nearest_transfer.get_position())) {
            std::cerr << "A nearest-momentum lookup returned a point outside its grid\n";
            return 1;
        }
        int closest_fermion_distance = std::numeric_limits<int>::max();
        int closest_transfer_distance = std::numeric_limits<int>::max();
        for (const auto& point : MomentumGrid::fermion_points()) {
            closest_fermion_distance =
                std::min(closest_fermion_distance, distance_squared(position, point.momentum.get_position()));
        }
        for (const auto& point : MomentumGrid::transfer_points()) {
            closest_transfer_distance =
                std::min(closest_transfer_distance, distance_squared(position, point.momentum.get_position()));
        }
        if (distance_squared(position, nearest_fermion.get_position()) != closest_fermion_distance ||
            distance_squared(position, nearest_transfer.get_position()) != closest_transfer_distance) {
            std::cerr << "A momentum lookup did not return a nearest available point\n";
            return 1;
        }
    }

    InteractionDataFrame frame;
    const momentum_iterator<L> K(1);
    const momentum_iterator<L> P(3);
    const momentum_iterator<L> Q(5);
    const auto resolved_k = MomentumGrid::nearest_fermion(K);
    const auto resolved_p = MomentumGrid::nearest_fermion(P);
    const auto resolved_q = MomentumGrid::nearest_transfer(Q);
    frame(K, P, Q) = 2.5;
    if (frame(resolved_k, resolved_p, resolved_q) != 2.5) {
        std::cerr << "Interaction lookup did not snap K, P, and Q to their respective grids\n";
        return 1;
    }

    for (const std::size_t index : InteractionDataFrame::active_independent_indices()) {
        const int k = static_cast<int>(index / (N * N));
        const int p = static_cast<int>((index / N) % N);
        const int q = static_cast<int>(index % N);
        if (!MomentumGrid::is_fermion(k) || !MomentumGrid::is_fermion(p) || !MomentumGrid::is_transfer(q)) {
            std::cerr << "An inactive interaction tuple was selected for flow evaluation\n";
            return 1;
        }
    }
    return 0;
}
