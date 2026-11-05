#pragma once

namespace NickelCUT {
#ifndef NICKELCUT_L
#define NICKELCUT_L 10
#endif

static constexpr int L = NICKELCUT_L;
static_assert(L >= 2 && L % 2 == 0, "NickelCUT lattice size L must be a positive even number.");
static constexpr int N = L * L;
}  // namespace NickelCUT
