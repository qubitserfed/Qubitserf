#include <algorithm>
#include <functional>
#include <vector>

#include "linear_algebra.hpp"
#include "combinatorics.hpp"

// input:  a number n and a function whose sole argument is a vector of bools
// effect: the function gets called on all length n bitstrings in lexicographical order
void partitions(int n, std::function<void(std::vector<bool>)> f) {
    std::vector<bool> stack;
    stack.reserve(n);

    std::function<void()> bkt = [&]() {
        if (stack.size() == n) {
            f(stack);
            return;
        }

        stack.push_back(0);
        bkt();
        stack.pop_back();

        stack.push_back(1);
        bkt();
        stack.pop_back();
    };

    bkt();
}

// input:  a number n and a function whose sole argument is a vector of bools
// effect: the function gets called on all length n bitstrings of Hamming weight k in lexicographical order
void combinations(int n, int k, std::function<void(std::vector<bool>) > f) {
    std::vector<bool> stack;
    stack.reserve(n);

    std::function<void(int)> bkt = [&](int weight) {
        if (stack.size() == n) {
            f(stack);
            return;
        }

        if (weight <= n - (int(stack.size()) + 1)) {
            stack.push_back(0);
            bkt(weight);
            stack.pop_back();
        }

        if (weight > 0) {
            stack.push_back(1);
            bkt(weight - 1);
            stack.pop_back();
        }
    };

    bkt(k);
}