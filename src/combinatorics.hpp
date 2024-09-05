#pragma once

#include <algorithm>
#include <functional>
#include <vector>

#include "linear_algebra.hpp"

void    partitions     (int, std::function<void(std::vector<bool>)>);
void    combinations   (int n, int k, std::function<void(std::vector<bool>) >);
