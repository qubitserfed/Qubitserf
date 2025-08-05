#include <algorithm>
#include <functional>
#include <iostream>
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

bool parallel_combinations(int n, int k, std::function<bool(BVector &)> f, int no_threads) {
    std::vector<u64> delimiters;
    std::vector<std::future<bool>> threads;

    u64 no_combinations = 1;
    for (int i = 1; i <= k; ++i)
        no_combinations = no_combinations * (n - i + 1) / i;

    for (int i = 0; i <= no_threads; ++i)
        delimiters.push_back(i * no_combinations / no_threads);

    ith_lexicographic_permutation(1, 1, 0); // initializes the internal static table, do not remove
    for (int i = 1; i < delimiters.size(); ++i) {
        int iv_start = delimiters[i - 1];
        int iv_end = delimiters[i];

        threads.push_back(std::async(std::launch::async, [n, k, iv_start, iv_end, f] () -> bool {
            BVector comb = ith_lexicographic_permutation(n, k, iv_start);
            for (u64 it = 0; it < iv_end - iv_start; ++it) {
                if (f(comb))
                    return true;
                advance_iterator(comb);
            }
            return false;
        }));
    }

    // todo: rather than wait for all threads to finish, we could return as soon as one thread finds a solution
    for (std::future<bool> &th: threads)
        if (th.get())
            return true;
    return false;
}

bool parallel_symplectic_combinations(int n, int k, std::function<bool(BVector &)> f, int no_threads) {
    return parallel_combinations(n, k, [&] (BVector &comb) {
        BVector symp_comb(2 * n);
        std::vector<int> pos;
        for (int i = 0; i < n; ++i)
            if (comb.get(i))
                pos.push_back(i);

        std::function<bool(int)> bkt = [&](int idx) -> bool {
            if (idx == pos.size())
                return f(symp_comb);

            symp_comb.set(2 * pos[idx], 1);
            symp_comb.set(2 * pos[idx] + 1, 0);
            if (bkt(idx + 1))
                return true;
    

            symp_comb.set(2 * pos[idx], 0);
            symp_comb.set(2 * pos[idx] + 1, 1);
            if (bkt(idx + 1))
                return true;

            symp_comb.set(2 * pos[idx], 1);
            symp_comb.set(2 * pos[idx] + 1, 1);
            if (bkt(idx + 1))
                return true;
            
            return false;
        };

        return bkt(0);
    }, no_threads);
}

void symplectic_combinations(int n, int k, std::function<void(std::vector<bool> &) > f) {
    std::vector<bool> stack;
    stack.reserve(2 * n);

    std::function<void(int)> bkt = [&](int weight) {
        if (stack.size() / 2 == n) {
            f(stack);
            return;
        }

        if (weight <= n - (int(stack.size()) / 2 + 1)) {
            stack.push_back(0);
            stack.push_back(0);
            bkt(weight);
            stack.pop_back();
            stack.pop_back();
        }

        if (weight > 0) {
            stack.push_back(1);
            stack.push_back(0);
            bkt(weight - 1);
            stack.pop_back();
            stack.pop_back();

            stack.push_back(0);
            stack.push_back(1);
            bkt(weight - 1);
            stack.pop_back();
            stack.pop_back();

            stack.push_back(1);
            stack.push_back(1);
            bkt(weight - 1);
            stack.pop_back();
            stack.pop_back();
        }
    };

    bkt(k);
}

BVector ith_lexicographic_permutation(int n, int k, u64 num) {
    const u64 infty = 0x3f3f3f3f3f3f3f3fULL;

    static bool initialized = false;
    static u64 c[200][200];

    if (!initialized) {
        for (int i = 0; i < 200; ++i)
            c[i][0] = 1;
        for (int i = 1; i < 200; ++i)
            for (int j = 1; j <= i; ++j)
                c[i][j] = std::min(infty, c[i - 1][j] + c[i - 1][j - 1]);
        initialized = true;
    }

    std::vector<bool> res;
    for (int i = 0; i < n; ++i) { // VERIFY
        if (c[n - i - 1][k] > num) {
            res.push_back(0);
        }
        else {
            num-= c[n - i - 1][k];
            k-= 1;
            res.push_back(1);
        }
    }
    my_assert(k == 0);

    return BVector(res);
}

bool advance_iterator(BVector &it) {
    const int n = it.n;
    auto &words = it.vec;

    if ((words[0] & -words[0]) + words[0] != 0) {
        using u64 = uint64_t;
        u64 x = words[0];
        if (x == 0)              // no 1's ⇒ no next
            return false;

        // mask off any bits above n
        u64 mask = (n < 64) ? ((u64(1) << n) - 1) : ~u64(0);
        x &= mask;

        u64 u = x & -x;          // lowest 1-bit
        u64 v = x + u;           // bump that bit
        if (v & ~mask)           // overflowed past bit n-1?
            return false;

        // spread the trailing 1's back down
        u64 w = v | (((v ^ x) / u) >> 2);
        w &= mask;

        words[0] = w;
        return true;
    }

    int bit = n, one_count = 0;
    while (bit > 0) {
        if (it.get(bit)) {
            one_count++;
        } else {
            if (one_count > 0) {
                bit++;
                break;
            }
        }
        bit--;
    }
    if (bit == 0)
        return false;

    // swap at (bit-1, bit)
    it.set(bit-1, true);
    it.set(bit,   false);

    // clear below bit and then set the last one_count-1 bits at the very end
    for (int i = bit; i < n; ++i)
        it.set(i, false);
    for (int i = n - 1; i >= n - one_count + 1; --i)
        it.set(i, true);

    return true;
}