#pragma once

#include <chrono>
#include <iostream>
#include <string>
#include <utility>
#include <memory>

#include "linear_algebra.hpp"
#include "combinatorics.hpp"


struct Timestamp {};
struct ResetTime {};

struct Printer {
    using clock = std::chrono::steady_clock;
    using time_point = std::chrono::steady_clock::time_point;

    time_point last_print_time;
    bool silent;

    void operator() () {
    }

    template<typename... Args>
    void operator()(int num, Args&&... args) {
        if (!silent) {
            std::cout << num;
        }
        ((*this)(std::forward<Args>(args)), ...);
    }

    template<typename... Args>
    void operator()(Timestamp, Args&&... args) {
        if (silent) return;
        auto now = clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_print_time);
        if (elapsed.count() >= 1000) {
            std::cout << "[" << (elapsed.count() / 1000.0) << "s] ";
        } else {
            std::cout << "[" << elapsed.count() << "ms] ";
        }
        last_print_time = now;
        ((*this)(std::forward<Args>(args)), ...);
    }

    template<typename... Args>
    void operator()(ResetTime, Args&&... args) {
        last_print_time = clock::now();
        ((*this)(std::forward<Args>(args)), ...);
    }

    template<typename... Args>
    void operator()(std::string s, Args&&... args) {
        if (!silent) {
            std::cout << s;
        }
        ((*this)(std::forward<Args>(args)), ...);
    }

    template<typename... Args>
    void operator()(double s, Args&&... args) {
        if (!silent) {
            std::cout << s;
        }
        ((*this)(std::forward<Args>(args)), ...);
    }

    Printer() : last_print_time(clock::now()), silent(false) {}
    Printer(bool not_silent) : last_print_time(clock::now()), silent(!not_silent) {}
};

using u32 = unsigned int;
using u64 = unsigned long long;

struct ParallelHashTable {
    std::vector<std::pair<BVector, BVector>> *table;
    std::unique_ptr<std::mutex[]> mutexes;

    int no_buckets, MASK;

    // hash function from https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
    u32 hashfn(u64);

    u32 hashfn(const BVector &);

    ParallelHashTable();

    bool insert(const BVector &key, const BVector &value);

    bool find(const BVector &key, const BVector &value);

    void reset(int pow2);

    void reset_for_code(int, int, bool);

    void clear();
};
