#pragma once

#include <algorithm>
#include <functional>
#include <future>
#include <thread>
#include <vector>

#include "linear_algebra.hpp"


struct COMPUTE_TYPE {
    bool CPU;
    bool GPU;
    int no_threads;
};

const COMPUTE_TYPE singlethreaded = {true, false, 1};
const COMPUTE_TYPE gpu = {false, true, 0};

void    partitions                      (int, std::function<void(std::vector<bool>)>);
void    combinations                    (int n, int k, std::function<void(std::vector<bool>) >);
void    symplectic_combinations         (int n, int k, std::function<void(std::vector<bool>&) >);
BVector ith_lexicographic_permutation   (int, int, u64);
bool    advance_iterator                (BVector &);
bool    parallel_combinations           (int n, int k, std::function<bool(BVector &)>, int no_threads = 1);
