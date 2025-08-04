#include <iostream>
#include <algorithm>
#include <map>
#include <mutex>
#include <set>
#include <utility>


#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "quantum_utilities.hpp"
#include "middle_algorithm.hpp"

using u32 = unsigned int;

struct SyndromeCell {
    int population_count;
    std::map< BVector, BVector > anticomm_syndromes;

    SyndromeCell() {
        population_count = 0;
    }
};

int css_middle_algorithm(BMatrix stab_mat, BMatrix anticomms) {
    const int m = stab_mat.m;

    std::map< BVector, std::set<BVector> > s0, s1;
    s0[BVector(stab_mat.n)].insert(BVector(anticomms.n));

    for (int d = 1; d <= m; ++d) {
        std::swap(s0, s1);
        s0.clear();

        try {
            combinations(m, d, [&] (const std::vector<bool> &v0) {
                BVector stab_syndome, anticomms_syndrome, v;
                std::map< BVector, std::set<BVector> >::iterator it;

                v = v0;
                stab_syndome = transposed_product(v, stab_mat);
                anticomms_syndrome = transposed_product(v, anticomms);

                // check for compatilble words of weight d-1
                it = s1.find(stab_syndome);
                if (it != s1.end()) {
                    if (it->second.size() >= 2) {
                        throw 2 * d - 1;
                    }
                    else if (*it->second.begin() != anticomms_syndrome) {
                        throw 2 * d - 1;
                    }
                }
            });

            combinations(m, d, [&] (const std::vector<bool> &v0) {
                BVector stab_syndome, anticomms_syndrome, v;
                std::map< BVector, std::set<BVector> >::iterator it;

                v = v0;
                stab_syndome = transposed_product(v, stab_mat);
                anticomms_syndrome = transposed_product(v, anticomms);

                it = s0.find(stab_syndome);
                if (it != s0.end()) {
                    if (it->second.size() >= 2) {
                        throw 2 * d;
                    }
                    else if (*it->second.begin() != anticomms_syndrome) {
                        throw 2 * d;
                    }
                    else {
                        it->second.insert(anticomms_syndrome);
                    }
                }
                else {
                    s0[stab_syndome].insert(anticomms_syndrome);
                }
            });
        }
        catch (int result) {
            return result;
        }
    }

    return 1;
}

int middle_algorithm(BMatrix stab_mat, BMatrix anticomms) {
    const int m = stab_mat.m / 2;

    std::vector<std::pair<u64, u64>> s0, s1;

    s0.emplace_back(0, 0);

    for (int d = 1; d <= m; ++d) {
        std::swap(s0, s1);
        s0.clear();
        long long reserve_size = 1;
        for (int i = 1; i <= d; ++i)
            reserve_size = reserve_size * (m - i + 1) / i;

        s0.reserve(reserve_size);

        try {
            long long cnt = 0;
            symplectic_combinations(m, d, [&] (std::vector<bool> &v0) {
                BVector stab_syndome, anticomms_syndrome, v;
                std::vector<std::pair<u64, u64>>::iterator it;

                v = v0;
                stab_syndome = BVector(stab_mat.n);
                for (int i = 0; i < stab_mat.n; ++i)
                    stab_syndome.set(i, sym_prod(stab_mat.row(i), v));

                anticomms_syndrome = BVector(anticomms.n);
                for (int i = 0; i < anticomms.n; ++i)
                    anticomms_syndrome.set(i, sym_prod(anticomms.row(i), v));

                // check for compatilble words of weight d-1
                it = lower_bound(s1.begin(), s1.end(), std::make_pair(stab_syndome.vec[0], 0ULL));
                if (it != s1.end() && it->first == stab_syndome.vec[0] && it->second != anticomms_syndrome.vec[0])
                    throw 2 * d - 1;

                if (it != s1.end())
                    it++;
                if (it != s1.end() && it->first == stab_syndome.vec[0] && it->second != anticomms_syndrome.vec[0])
                    throw 2 * d - 1;
                cnt++;
            });

            symplectic_combinations(m, d, [&] (const std::vector<bool> &v0) {
                BVector stab_syndome, anticomms_syndrome, v;
                std::vector<std::pair<u64, u64>>::iterator it;

                v = v0;
                stab_syndome = BVector(stab_mat.n);
                for (int i = 0; i < stab_mat.n; ++i)
                    stab_syndome.set(i, sym_prod(stab_mat.row(i), v));

                anticomms_syndrome = BVector(anticomms.n);
                for (int i = 0; i < anticomms.n; ++i)
                    anticomms_syndrome.set(i, sym_prod(anticomms.row(i), v));

                s0.emplace_back(stab_syndome.vec[0], anticomms_syndrome.vec[0]);
            });

            std::sort(s0.begin(), s0.end());

            symplectic_combinations(m, d, [&] (const std::vector<bool> &v0) {
                BVector stab_syndome, anticomms_syndrome, v;
                std::vector<std::pair<u64, u64>>::iterator it;

                v = v0;
                stab_syndome = BVector(stab_mat.n);
                for (int i = 0; i < stab_mat.n; ++i)
                    stab_syndome.set(i, sym_prod(stab_mat.row(i), v));

                anticomms_syndrome = BVector(anticomms.n);
                for (int i = 0; i < anticomms.n; ++i)
                    anticomms_syndrome.set(i, sym_prod(anticomms.row(i), v));

                // check for compatilble words of weight d-1
                it = lower_bound(s0.begin(), s0.end(), std::make_pair(stab_syndome.vec[0], 0ULL));
                if (it != s0.end() && it->first == stab_syndome.vec[0] && it->second != anticomms_syndrome.vec[0])
                    throw 2 * d;

                it++;
                if (it != s0.end() && it->first == stab_syndome.vec[0] && it->second != anticomms_syndrome.vec[0])
                    throw 2 * d;
            });
        }
        catch (int result) {
            return result;
        }
    }

    return 1;
}

std::pair<int, int> get_zx_distances_with_middle(BMatrix stab_mat, bool verbose_flag) {
    BMatrix closure_mat, x_stab, z_stab, x_ops, z_ops;

    closure_mat = logical_operators(stab_mat);

    std::tie(z_ops, x_ops) = zx_parts(closure_mat);
    std::tie(z_stab, x_stab) = zx_parts(stab_mat);

    to_row_echelon(z_stab);
    to_row_echelon(x_stab);

    z_stab.remove_zeros();
    x_stab.remove_zeros();

    z_ops.remove_zeros();
    x_ops.remove_zeros();

    const int z_dist = css_middle_algorithm(x_stab, x_ops);
    const int x_dist = css_middle_algorithm(z_stab, z_ops);

    return std::make_pair(z_dist, x_dist);
}


const long long MAX_BUCKETS = 1LL << 25;

u32 hashfn(u64 key) {
    key = key * 0xbf58476d1ce4e5b9ULL;
    key = key ^ (key >> 32);
    key = key * 0xbf58476d1ce4e5b9ULL;
    key = key ^ (key >> 32);
    key = key * 0xbf58476d1ce4e5b9ULL;
    if (MAX_BUCKETS == 1LL << 32) // compile time constant
        return key;
    else
        return key % MAX_BUCKETS;
}

u32 hashfn(const BVector &v) {
    u64 key = 0;
    for (int i = 0; i < v.vec.size(); ++i)
        key = key + hashfn(v.vec[i]);
    return key;
}

struct ParallelHashTable {
    std::vector<std::pair<BVector, BVector>> table[MAX_BUCKETS];
    std::vector<std::mutex> mutexes;
    
    ParallelHashTable() : mutexes(MAX_BUCKETS) {}

    bool insert(const BVector &key, const BVector &value) { // returns true if the key was already present with a different value
        u32 hash = hashfn(key) % MAX_BUCKETS;
        std::lock_guard<std::mutex> lock(mutexes[hash]);
        for (auto &pair : table[hash]) {
            if (pair.first == key) {
                if (pair.second != value) {
                    return true;
                }
                return false;
            }
        }
        table[hash].emplace_back(key, value);
        return false;
    }

    bool find(const BVector &key, const BVector &value) {
        u32 hash = hashfn(key) % MAX_BUCKETS;
        std::lock_guard<std::mutex> lock(mutexes[hash]);
        for (auto &pair : table[hash]) {
            if (pair.first == key && pair.second != value) {
                return true;
            }
        }
        return false;
    }

    void clear() {
        for (int i = 0; i < MAX_BUCKETS; ++i)
            table[i].clear();
    }
};

int parallel_middle_algorithm(BMatrix stab_mat, BMatrix anticomms, COMPUTE_TYPE compute_type) {
    const int m = stab_mat.m / 2;

    ParallelHashTable *s0, *s1;
    s0 = new ParallelHashTable();
    s1 = new ParallelHashTable();

    s0->insert(BVector(stab_mat.n), BVector(anticomms.n));

    for (int d = 1; d <= m; ++d) {
        std::swap(s0, s1);
        s0->clear();

        bool res = false;

        res = parallel_symplectic_combinations(m, d, [&] (BVector &v) -> bool {
            BVector stab_syndome, anticomms_syndrome;

            stab_syndome = BVector(stab_mat.n);
            for (int i = 0; i < stab_mat.n; ++i)
                stab_syndome.set(i, sym_prod(stab_mat.row(i), v));

            anticomms_syndrome = BVector(anticomms.n);
            for (int i = 0; i < anticomms.n; ++i)
                anticomms_syndrome.set(i, sym_prod(anticomms.row(i), v));

            // check for compatilble words of weight d-1
            if (s1->find(stab_syndome, anticomms_syndrome))
                return true;
            s0->insert(stab_syndome, anticomms_syndrome);
            return false;
        }, compute_type.no_threads);

        if (res) {
            delete s0;
            delete s1;
            return 2 * d - 1;
        }

        res = parallel_symplectic_combinations(m, d, [&] (BVector &v) -> bool {
            BVector stab_syndome, anticomms_syndrome;

            stab_syndome = BVector(stab_mat.n);
            for (int i = 0; i < stab_mat.n; ++i)
                stab_syndome.set(i, sym_prod(stab_mat.row(i), v));

            anticomms_syndrome = BVector(anticomms.n);
            for (int i = 0; i < anticomms.n; ++i)
                anticomms_syndrome.set(i, sym_prod(anticomms.row(i), v));

            if (s0->find(stab_syndome, anticomms_syndrome))
                return true;
            return false;
        }, compute_type.no_threads);

        if (res) {
            delete s0;
            delete s1;
            return 2 * d;
        }
    }

    delete s0;
    delete s1;
    return 1;
}

int parallel_css_middle_algorithm(BMatrix stab_mat, BMatrix anticomms, COMPUTE_TYPE compute_type) {
    const int m = stab_mat.m;

    ParallelHashTable *s0, *s1;
    s0 = new ParallelHashTable();
    s1 = new ParallelHashTable();
    s0->insert(BVector(stab_mat.n), BVector(anticomms.n));

    for (int d = 1; d <= m; ++d) {
        std::swap(s0, s1);
        s0->clear();

        bool found = false;
        found = parallel_combinations(m, d, [&] (BVector &v) {
            BVector stab_syndome = transposed_product(v, stab_mat);
            BVector anticomms_syndrome = transposed_product(v, anticomms);

            if (s1->find(stab_syndome, anticomms_syndrome)) {
                return true;
            }
            return false;
        }, compute_type.no_threads);

        if (found) {
            delete s0;
            delete s1;
            return 2 * d - 1;
        }

        found = parallel_combinations(m, d, [&] (BVector &v) {
            BVector stab_syndome = transposed_product(v, stab_mat);
            BVector anticomms_syndrome = transposed_product(v, anticomms);
            return s0->insert(stab_syndome, anticomms_syndrome);
        }, compute_type.no_threads);

        if (found) {
            delete s0;
            delete s1;
            return 2 * d;
        }
    }

    delete s0;
    delete s1;
    return 1;
}

std::pair<int, int> get_zx_distances_with_parallelized_middle(BMatrix stab_mat, COMPUTE_TYPE compute_type) {
    BMatrix closure_mat, x_stab, z_stab, x_ops, z_ops;

    closure_mat = logical_operators(stab_mat);

    std::tie(z_ops, x_ops) = zx_parts(closure_mat);
    std::tie(z_stab, x_stab) = zx_parts(stab_mat);

    to_row_echelon(z_stab);
    to_row_echelon(x_stab);

    z_stab.remove_zeros();
    x_stab.remove_zeros();

    z_ops.remove_zeros();
    x_ops.remove_zeros();

    const int z_dist = parallel_css_middle_algorithm(x_stab, x_ops, compute_type);
    const int x_dist = parallel_css_middle_algorithm(z_stab, z_ops, compute_type);

    return std::make_pair(z_dist, x_dist);
}

int get_distance_with_parallelized_middle(BMatrix stab_mat, COMPUTE_TYPE compute_type) {
    if (is_css(stab_mat)) {
        int z_dist, x_dist;
        std::tie(z_dist, x_dist) = get_zx_distances_with_parallelized_middle(stab_mat, compute_type);
        return std::min(z_dist, x_dist);
    }
    else {
        return parallel_middle_algorithm(stab_mat, logical_operators(stab_mat), compute_type);
    }
}

int get_distance_with_middle(BMatrix stab_mat) {
    if (is_css(stab_mat)) {
        int z_dist, x_dist;
        std::tie(z_dist, x_dist) = get_zx_distances_with_middle(stab_mat);
        return std::min(z_dist, x_dist);
    }
    else {
        return middle_algorithm(stab_mat, logical_operators(stab_mat));
    }
}
