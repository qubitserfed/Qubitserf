#pragma once

#include <algorithm>
#include <utility>

#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "quantum_utilities.hpp"
#include "utility.hpp"

enum ALGORITHMS {
    MIDDLE_ALGORITHM,
    BROUWER_ZIMMERMAN_ALGORITHM,
};

std::pair<int, int> get_zx_distances(BMatrix, ALGORITHMS, COMPUTE_TYPE, bool = false);
int get_distance(BMatrix, ALGORITHMS, COMPUTE_TYPE, bool = false);
