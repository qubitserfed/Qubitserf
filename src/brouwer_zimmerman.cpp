#include <iostream>
#include <algorithm>
#include <vector>

#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "quantum_utilities.hpp"
#include "brouwer_zimmerman.hpp"

// input:  two *sorted* int vectors
// output: a sorted vector containing the elements of a that are not contained in b
std::vector<int> list_diff(std::vector<int> a, std::vector<int> b) {
    std::vector<int> res;
    int i = 0;
    int j = 0;

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

// input: generator matrix of an arbitrary code
// output: distance of the code
int brouwer_zimmerman(BMatrix stab_mat, BMatrix code_mat) {
    const int n = code_mat.n;
    const int m = code_mat.m;

    std::vector< std::pair<BMatrix, int> > gamma_seq = brouwer_zimmerman_sequence(code_mat);
    int inner_bound = 2e9, outer_bound = -2e9;

//    std::cout << "BZS:\n" << std::endl;
//    for (auto mat: gamma_seq)
//        print(mat.first);

    for (int d = 1; d <= n; ++d) {
        for (auto gen_pair: gamma_seq) {
            BMatrix gamma_mat = gen_pair.first;
            BMatrix transposed_gamma_mat = transpose(gamma_mat);

            combinations(n, d, [&](std::vector<bool> v0) {
                BVector vec(v0);

                BVector codeword = transposed_product(vec, transposed_gamma_mat);

                if (!in_span(stab_mat, codeword))
                    inner_bound = std::min(inner_bound, codeword.weight());
            });
        }

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


std::pair<int, int> get_zx_distances(BMatrix stab_mat) {
    BMatrix z_closed_mat, x_closed_mat, x_stab, z_stab, x_closed, z_closed, junk;

    to_row_echelon(stab_mat);

//    print(stab_mat);
    
    z_closed_mat = prefered_isotropic_closure(stab_mat, false);
    x_closed_mat = prefered_isotropic_closure(stab_mat, true);
    
    std::tie(z_stab, x_stab) = zx_parts(stab_mat);
//    print(z_stab);

    x_stab.remove_zeros();
    z_stab.remove_zeros();
    
    to_row_echelon(x_stab);
    to_row_echelon(z_stab);

    std::tie(z_closed, junk) = zx_parts(z_closed_mat);
    std::tie(junk, x_closed) = zx_parts(x_closed_mat);

    z_stab.remove_zeros(); z_stab.sort_rows();
    z_closed.remove_zeros(); z_closed.sort_rows();

    x_stab.remove_zeros(); x_stab.sort_rows();
    x_closed.remove_zeros(); x_closed.sort_rows();

//    print(z_stab);
//    print(z_closed);

    const int z_dist = brouwer_zimmerman(z_stab, z_closed);
    const int x_dist = brouwer_zimmerman(x_stab, x_closed);

    return std::make_pair(z_dist, x_dist);
}


int get_distance(BMatrix stab_mat) {
    int z_dist, x_dist;
    std::tie(z_dist, x_dist) = get_zx_distances(stab_mat);
    return std::min(z_dist, x_dist);
}
