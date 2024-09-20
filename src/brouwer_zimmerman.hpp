#pragma once

#include <algorithm>
#include <vector>

#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "quantum_utilities.hpp"

enum COMPUTE_TYPE {
    CPU,
    GPU
};

std::vector<int>                        list_diff                   (std::vector<int>, std::vector<int>);
std::vector<int>                        zero_columns                (BMatrix, std::vector<int>);
std::vector<std::pair<BMatrix, int> >   brouwer_zimmerman_sequence  (BMatrix);
int                                     brouwer_zimmerman           (BMatrix, BMatrix);
std::pair<int, int>                     get_zx_distances            (BMatrix, COMPUTE_TYPE = CPU);
int                                     get_distance                (BMatrix, COMPUTE_TYPE = CPU);
