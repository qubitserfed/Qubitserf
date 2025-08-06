#pragma once

#include <algorithm>
#include <vector>

#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "quantum_utilities.hpp"

std::vector<std::pair<BMatrix, int> >   brouwer_zimmerman_sequence  (BMatrix);
int                                     brouwer_zimmerman           (BMatrix, BMatrix, COMPUTE_TYPE = singlethreaded, bool = false);
std::pair<int, int>                     get_zx_distances            (BMatrix, COMPUTE_TYPE = singlethreaded, bool = false);
