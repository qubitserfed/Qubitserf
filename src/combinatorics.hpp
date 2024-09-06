#pragma once

#include <algorithm>
#include <functional>
#include <future>
#include <thread>
#include <vector>

#include "linear_algebra.hpp"

void    partitions     (int, std::function<void(std::vector<bool>)>);
void    combinations   (int n, int k, std::function<void(std::vector<bool>) >);

int parallel_combinations(
    int,
    int,
    std::function< int(std::vector<bool>) >,
    std::function< int(int, int) >,
    int
);
