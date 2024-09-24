#include <algorithm>
#include <future>
#include <iostream>
#include <utility>
#include <vector>

#include "parallel_bz.hpp"

int exponential_part_parallel(BMatrix stab_mat, std::vector<BMatrix> transposed_gamma_sequence, int n, int d, int no_threads) {
    std::vector<std::future<int>> threads;
    std::vector<u64> delimiters;
    u64 no_combinations = 1;

    for (int i = 1; i <= d; ++i)
        no_combinations = no_combinations * (n - i + 1) / i;

    for (int i = 0; i <= no_threads; ++i)
        delimiters.push_back(i * no_combinations / no_threads);

    ith_lexicographic_permutation(1, 1, 0);
    for (int i = 1; i < delimiters.size(); ++i) {
        int iv_start = delimiters[i - 1];
        int iv_end = delimiters[i];

        threads.push_back(std::async([transposed_gamma_sequence, stab_mat, n, d, iv_start, iv_end] () -> int {
            BVector comb = ith_lexicographic_permutation(n, d, iv_start);
            int bound  = 2e9;
            for (u64 it = 0; it < iv_end - iv_start; ++it) {
                for (auto transposed_gamma_mat: transposed_gamma_sequence) {
                    BVector codeword = transposed_product(comb, transposed_gamma_mat);
                    if (!in_span(stab_mat, codeword))
                        bound = std::min(bound, codeword.weight());
                }
                advance_iterator(comb);
            }
            return bound;
        }));
    }

    int bound = 2e9;
    for (std::future<int> &th: threads)
        bound = std::min(bound, th.get());

    return bound;
}
