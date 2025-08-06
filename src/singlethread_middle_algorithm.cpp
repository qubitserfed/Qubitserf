#include <iostream>
#include <algorithm>
#include <map>
#include <mutex>
#include <set>
#include <utility>

#include "utility.hpp"
#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "quantum_utilities.hpp"
#include "singlethread_middle_algorithm.hpp"

using u32 = unsigned int;

struct SyndromeCell {
    int population_count;
    std::map< BVector, BVector > anticomm_syndromes;

    SyndromeCell() {
        population_count = 0;
    }
};

int css_middle_algorithm(BMatrix stab_mat, BMatrix anticomms, bool verbose_flag) {
    const int m = stab_mat.m;
    Printer printer(verbose_flag);

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

            printer("Distance bound: >", 2 * d - 1, "\nElapsed:", Timestamp(), "\n", ResetTime());

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

            printer("Distance bound: >", 2 * d, "\nElapsed:", Timestamp(), "\n", ResetTime());
        }
        catch (int result) {
            printer("Distance bound: =", result, "\nElapsed:", Timestamp(), "\n", ResetTime());
            return result;
        }

    }

    return 1;
}

int middle_algorithm(BMatrix stab_mat, BMatrix anticomms, bool verbose_flag) {
    const int m = stab_mat.m / 2;

    Printer printer(verbose_flag);
    std::vector<std::pair<u64, u64>> s0, s1;

    s0.emplace_back(0, 0);
    printer(ResetTime());

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

            printer("Distance bound: >", 2 * d - 1, "\nElapsed:", Timestamp(), "\n", ResetTime());

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

                if (it != s0.end())
                    it++;
                if (it != s0.end() && it->first == stab_syndome.vec[0] && it->second != anticomms_syndrome.vec[0])
                    throw 2 * d;
            });

            printer("Distance bound: >", 2 * d, "\nElapsed:", Timestamp(), "\n", ResetTime());
        }
        catch (int result) {
            printer("Distance: =", result, "\nElapsed:", Timestamp(), "\n", ResetTime());
            return result;
        }

    }

    return 1;
}

std::pair<int, int> get_zx_distances_with_middle(BMatrix stab_mat, bool verbose_flag) {
    BMatrix closure_mat, x_stab, z_stab, x_ops, z_ops;
    Printer printer(verbose_flag);

    closure_mat = logical_operators(stab_mat);

    std::tie(z_ops, x_ops) = zx_parts(closure_mat);
    std::tie(z_stab, x_stab) = zx_parts(stab_mat);

    to_row_echelon(z_stab);
    to_row_echelon(x_stab);

    z_stab.remove_zeros();
    x_stab.remove_zeros();

    z_ops.remove_zeros();
    x_ops.remove_zeros();

    printer("Z-distance:\n");
    const int z_dist = css_middle_algorithm(x_stab, x_ops, verbose_flag);

    printer("X-distance:\n");
    const int x_dist = css_middle_algorithm(z_stab, z_ops, verbose_flag);

    return std::make_pair(z_dist, x_dist);
}

int get_distance_with_middle(BMatrix stab_mat, bool verbose_flag) {
    if (is_css(stab_mat)) {
        int z_dist, x_dist;
        std::tie(z_dist, x_dist) = get_zx_distances_with_middle(stab_mat, verbose_flag);
        return std::min(z_dist, x_dist);
    }
    else {
        return middle_algorithm(stab_mat, logical_operators(stab_mat), verbose_flag);
    }
}

