#include <algorithm>
#include <functional>
#include <vector>

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

int parallel_combinations(
    int n,
    int k,
    std::function< int(std::vector<bool>) > f,
    std::function< int(int,int) > merge,
    int NO_THREADS
) {
    int c[64][64]; // c[i][j] = min(1024, i choose j)
    std::vector<bool> stack;
    stack.reserve(n);

    // compute c[i][j] using Pascal's triangle
    for (int i = 0; i < 64; ++i)
        c[i][0] = 1;
    for (int i = 1; i < 64; ++i) {
        for (int j = 1; j < 64; ++j) {
            c[i][j] = std::min(1024, c[i - 1][j] + c[i - 1][j - 1]);
        }
    }

    std::function<int(int, std::vector<bool>)> parallel_bkt = [&](int weight, std::vector<bool> stack) -> int {
        if (stack.size() == n) {
            return f(stack);
        }

        bool one_set = false;
        bool zero_set = false;
        std::future<int> one_val, zero_val;
        if (weight <= n - (int(stack.size()) + 1)) {
            stack.push_back(0);
    
            zero_val = std::async(std::launch::async, parallel_bkt, weight, stack);
    
            zero_set = true;
            stack.pop_back();
        }

        if (weight > 0) {
            stack.push_back(1);
            
            one_val = std::async(std::launch::async, parallel_bkt, weight - 1, stack);

            one_set = true;
            stack.pop_back();
        }

        if (one_set && zero_set)
            return merge(zero_val.get(), one_val.get());
        else if (one_set)
            return one_val.get();
        else if (zero_set)
            return zero_val.get();
        else {
            my_assert(false);
            return zero_val.get();
        }
    };

    std::function<int(int)> bkt = [&](int weight) -> int {
        if (stack.size() == n) {
            return f(stack);
        }

        bool one_set = false;
        bool zero_set = false;
        int one_val, zero_val;

        if (weight <= n - (int(stack.size()) + 1)) {
            stack.push_back(0);
    
            if (n - int(stack.size()) < 64 && c[n - int(stack.size())][weight] <= NO_THREADS)
                zero_val = parallel_bkt(weight, stack);
            else
                zero_val = bkt(weight);
    
            zero_set = true;
            stack.pop_back();
        }

        if (weight > 0) {
            stack.push_back(1);

            if (n - int(stack.size()) < 64 && c[n - int(stack.size())][weight - 1] <= NO_THREADS)
                one_val = parallel_bkt(weight, stack);
            else
                one_val = bkt(weight - 1);

            one_set = true;
            stack.pop_back();
        }

        if (one_set && zero_set)
            return merge(zero_val, one_val);
        else if (one_set)
            return one_val;
        else if (zero_set)
            return zero_val;
        else {
            return zero_val;
            my_assert(false);
        }
    };

    return bkt(k);
}
