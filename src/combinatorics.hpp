#pragma once

#include <algorithm>
#include <functional>
#include <future>
#include <thread>
#include <vector>

#include "linear_algebra.hpp"

void    partitions                      (int, std::function<void(std::vector<bool>)>);
void    combinations                    (int n, int k, std::function<void(std::vector<bool>) >);
void    symplectic_combinations         (int n, int k, std::function<void(std::vector<bool>&) >);
BVector ith_lexicographic_permutation   (int, int, u64);
bool    advance_iterator                (BVector &);
