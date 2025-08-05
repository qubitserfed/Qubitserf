#include <utility>
#include "linear_algebra.hpp"
#include "combinatorics.hpp"

int css_middle_algorithm(BMatrix, BMatrix, bool = false);
int middle_algorithm(BMatrix, BMatrix, bool = false);
int get_distance_with_middle(BMatrix, bool = false);
std::pair<int, int> get_zx_distances_with_middle(BMatrix, bool = false);
int parallel_css_middle_algorithm(BMatrix, BMatrix, COMPUTE_TYPE, bool = false);
int get_distance_with_parallelized_middle(BMatrix, COMPUTE_TYPE, bool = false);
std::pair<int, int> get_zx_distances_with_parallelized_middle(BMatrix, COMPUTE_TYPE, bool = false);
