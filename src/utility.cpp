#include "utility.hpp"

u32 ParallelHashTable::hashfn(u64 key) {
    key = key * 0xbf58476d1ce4e5b9ULL;
    key = key ^ (key >> 32);
    key = key * 0xbf58476d1ce4e5b9ULL;
    return key & MASK;
}

u32 ParallelHashTable::hashfn(const BVector &v) {
    u64 key = 0;
    for (int i = 0; i < v.vec.size(); ++i)
        key = hashfn(key + v.vec[i]);
    return key;
}

bool ParallelHashTable::insert(const BVector &key, const BVector &value) {
    u32 hash = hashfn(key);
    std::lock_guard<std::mutex> lock(mutexes[hash]);
    for (auto &pair : table[hash]) {
        if (pair.first == key) {
            if (pair.second != value) {
                return true;
            }
        }
    }
    table[hash].emplace_back(key, value);
    return false;
}

bool ParallelHashTable::find(const BVector &key, const BVector &value) {
    u32 hash = hashfn(key);
    std::lock_guard<std::mutex> lock(mutexes[hash]);
    for (auto &pair : table[hash]) {
        if (pair.first == key && pair.second != value) {
            return true;
        }
    }
    return false;
}

void ParallelHashTable::reset(int pow2) {
    if (table != nullptr) {
        delete[] table;
        table = nullptr;
    }
    no_buckets = 1 << pow2;
    MASK = no_buckets - 1;
    table = new std::vector<std::pair<BVector, BVector>>[no_buckets];
    mutexes = std::make_unique<std::mutex[]>(no_buckets);
}

void ParallelHashTable::reset_for_code(int n, int td, bool css) {
    double approx_combinations = 1;
    int bucket_pow2 = 0;
    for (int i = 1; i <= td; ++i)
        approx_combinations = approx_combinations * (n - i + 1) / i * 3;

    if (approx_combinations > 2e9) {
        bucket_pow2 = 25;
    }
    else {
        while (approx_combinations > 0.5) {
            approx_combinations /= 2; 
            bucket_pow2+= 1;
        }
        bucket_pow2 = std::max(16, std::min(25, bucket_pow2));
    }
    reset(bucket_pow2);
}

void ParallelHashTable::clear() {
    for (int i = 0; i < no_buckets; ++i)
        table[i].clear();
}

ParallelHashTable::ParallelHashTable() : table(nullptr) {
    reset(16);
}
