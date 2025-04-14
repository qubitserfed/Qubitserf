#include <utility>
#include "linear_algebra.hpp"
#include "combinatorics.hpp"

int css_middle_algorithm(BMatrix, BMatrix);
int middle_algorithm(BMatrix, BMatrix);
int get_distance_with_middle(BMatrix);
std::pair<int, int> get_zx_distances_with_middle(BMatrix);
int parallel_css_middle_algorithm(BMatrix, BMatrix, COMPUTE_TYPE);
int get_distance_with_parallelized_middle(BMatrix, COMPUTE_TYPE);
std::pair<int, int> get_zx_distances_with_parallelized_middle(BMatrix, COMPUTE_TYPE);
