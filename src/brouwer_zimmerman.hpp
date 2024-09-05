#pragma once

#include <algorithm>
#include <vector>

#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "quantum_utilities.hpp"

// input:  two *sorted* int vectors
// output: a sorted vector containing the elements of a that are not contained in b
std::vector<int>                        list_diff                   (std::vector<int>, std::vector<int>);
std::vector<int>                        zero_columns                (BMatrix, std::vector<int>);
std::vector<std::pair<BMatrix, int> >   brouwer_zimmerman_sequence  (BMatrix);
int                                     brouwer_zimmerman           (BMatrix, BMatrix);
std::pair<int, int>                     get_zx_distances            (BMatrix);
int                                     get_distance                (BMatrix);
