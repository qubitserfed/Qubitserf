#pragma once

#include "linear_algebra.hpp"
#include "combinatorics.hpp"

int                     get_operator_weight         (BMatrix, BVector, COMPUTE_TYPE, bool);
std::pair<int, int>     get_zx_operator_weight      (BMatrix, BVector, COMPUTE_TYPE, bool);
