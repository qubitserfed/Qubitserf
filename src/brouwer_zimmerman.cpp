#include <algorithm>
#include <future>
#include <iostream>
#include <thread>
#include <vector>

#include <cstdlib>


#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "quantum_utilities.hpp"
#include "brouwer_zimmerman.hpp"
#include "parallel_bz.hpp"

// input:  two int vectors
// output: a sorted vector containing the elements of a that are not contained in b
std::vector<int> list_diff(std::vector<int> a, std::vector<int> b) {
    std::vector<int> res;
    int i = 0;
    int j = 0;

    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());

    res.reserve(a.size());
    for (; i < a.size(); ++i) {
        while (j < b.size() && b[j] < a[i])
            ++j;
        if (j == b.size() || a[i] != b[j])
            res.push_back(a[i]);
    }

    return res;
}

// input:  a matrix M and a list of columns L
// output: the subset of columns in L which are zero
std::vector<int> zero_columns(BMatrix mat, std::vector<int> cols) {
    const int n = mat.n;
    const int m = mat.m;

    std::vector<int> result;
    for (int j: cols) {
        bool all_zero = true;
        
        for (int i = 0; all_zero && i < n; ++i)
            if (mat.get(i, j))
                all_zero = false;

        if (all_zero)
            result.push_back(j);
    }

    return result;
}

// input:  a generator matrix gen
// output: the sequence of matrices Gamma_0, Gamma_1, \dots, Gamma_n paired with the numbers r_0, r_1, \dots, r_n
std::vector<std::pair<BMatrix, int>> brouwer_zimmerman_sequence(BMatrix gen) {
    const int n = gen.n;
    const int m = gen.m;
    
    std::vector<std::pair<BMatrix, int>> gamma_seq;
    std::vector<int> active_columns;
    BMatrix working_matrix;
    
    active_columns.resize(m);
    for (int i = 0; i < m; ++i)
        active_columns[i] = i;
    
    working_matrix = gen;
    while (!active_columns.empty()) {
        std::vector<int> eliminated_columns = restricted_row_echelon(working_matrix, active_columns);
        std::vector<int> new_zero_columns = zero_columns(working_matrix, active_columns);

        gamma_seq.emplace_back(working_matrix, eliminated_columns.size());

        active_columns = list_diff(active_columns, eliminated_columns);
        active_columns = list_diff(active_columns, new_zero_columns);
    }

    return gamma_seq;
}

int exponential_part_cpu(BMatrix stab_mat, std::vector<BMatrix> transposed_gamma_sequence, int n, int d) {
    int bound = 2e9;

    for (auto transposed_gamma_mat: transposed_gamma_sequence) {
        combinations(
            n,
            d,
            [&](std::vector<bool> v0) -> void {
                BVector vec(v0);
                BVector codeword = transposed_product(vec, transposed_gamma_mat);

                if (!in_span(stab_mat, codeword))
                    bound = std::min(bound, codeword.weight());
            }
        );
    }

    return bound;
}

int brouwer_zimmerman(BMatrix stab_mat, BMatrix code_mat, COMPUTE_TYPE compute_type) {
    const int n = code_mat.n;
    const int m = code_mat.m;

    std::vector< std::pair<BMatrix, int> > gamma_seq = brouwer_zimmerman_sequence(code_mat);
    std::vector< BMatrix > transposed_gamma_seq;

    for (auto &entry_pair : gamma_seq)
        transposed_gamma_seq.push_back(transpose(entry_pair.first));

    int inner_bound = 2e9, outer_bound = -2e9;

    for (int d = 1; d <= n; ++d) {
        int weight_d_bound;

        if (compute_type.CPU) {
            weight_d_bound = compute_type.no_threads == 1 ?
                exponential_part_cpu(stab_mat, transposed_gamma_seq, n, d) :
                exponential_part_parallel(stab_mat, transposed_gamma_seq, n, d, compute_type.no_threads);
        }
        else if (compute_type.GPU) {
            my_assert(0);
        }
        else {
            my_assert(0);
        }

        inner_bound = std::min(inner_bound, weight_d_bound);

        outer_bound = 0;
        for (auto gen_pair: gamma_seq) {
            const int r = gen_pair.second;

            outer_bound+= std::max(0, (d + 1) - (n - r));
        }

        if (inner_bound <= outer_bound)
            break;
    }

    return inner_bound;

}

std::pair<int, int> get_zx_distances(BMatrix stab_mat, COMPUTE_TYPE compute_type) {
    BMatrix closure_mat, x_stab, z_stab, x_closed, z_closed;

    closure_mat = isotropic_closure(stab_mat);

    std::tie(z_closed, x_closed) = zx_parts(closure_mat);
    std::tie(z_stab, x_stab) = zx_parts(stab_mat);

    to_row_echelon(z_stab);
    to_row_echelon(x_stab);

    z_stab.remove_zeros();
    x_stab.remove_zeros();

    z_closed.remove_zeros();
    x_closed.remove_zeros();;

    const int z_dist = brouwer_zimmerman(z_stab, z_closed, compute_type);
    const int x_dist = brouwer_zimmerman(x_stab, x_closed, compute_type);

    return std::make_pair(z_dist, x_dist);
}


// input: generator matrix of an arbitrary code
// output: distance of the code
int get_distance(BMatrix stab_mat, COMPUTE_TYPE compute_type) {
    int z_dist, x_dist;
    std::tie(z_dist, x_dist) = get_zx_distances(stab_mat, compute_type);
    return std::min(z_dist, x_dist);
}
