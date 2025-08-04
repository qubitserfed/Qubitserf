#include <algorithm>
#include <fstream>
#include <iostream>
#include <functional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <cctype>

#include "linear_algebra.hpp"
#include "brouwer_zimmerman.hpp"
#include "middle_algorithm.hpp"

BMatrix steane_code() {
    const int n = 6;
    const int m = 7;

    // 0 is I, 1 is Z, 2 is X, 3 is Y
    const int entries[n][m] = {
        {0, 0, 0, 1, 1, 1, 1},
        {0, 0, 0, 2, 2, 2, 2},
        {0, 1, 1, 0, 0, 1, 1},
        {0, 2, 2, 0, 0, 2, 2},
        {1, 0, 1, 0, 1, 0, 1},
        {2, 0, 2, 0, 2, 0, 2}
    };

    BMatrix steane(n, 2 * m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (entries[i][j] & 1)
                steane.set(i, 2 * j, 1);
            if (entries[i][j] & 2)
                steane.set(i, 2 * j + 1, 1);
        }
    }

    return steane;
}

BMatrix bmatrix_conversion(std::vector<std::string> code) {
    const int m = 2 * code.back().size();
    const int n = code.size();

    BMatrix result;

    for (int i = 0; i < n; ++i) {
        BVector row(m);

        for (int j = 0; j < code[i].size(); ++j) {
            switch (code[i][j]) {
                case '_':
                case '.':
                case 'I':
                    break;
                case 'X':
                    row.set(2 * j + 1, 1);
                    break;
                case 'Z':
                    row.set(2 * j, 1);
                    break;
                case 'Y':
                    row.set(2 * j, 1);
                    row.set(2 * j + 1, 1);
                    break;
                default:
                    my_assert(0);
                    break;
            }
        }

        result.append_row(row);
    }

    return result;
}

int get_middle_distance(BMatrix stab_mat) {
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

    const int z_dist = middle_algorithm(x_stab, x_ops);
    const int x_dist = middle_algorithm(z_stab, z_ops);

    return std::min(z_dist, x_dist);
}

std::vector<std::tuple<BMatrix, int, int>> code_list(std::string filename) {
    std::ifstream fi(filename);

    std::vector<std::tuple<BMatrix, int, int>> bencodes;
    std::vector<std::string> code;
    std::string current;
    int low_bound, high_bound;

    while (getline(fi, current)) {
        if (current == "") {
            bencodes.push_back(std::make_tuple(bmatrix_conversion(code), low_bound, high_bound));
            low_bound = high_bound = 0;
            code.clear();
        }
        else {
            if (std::any_of(current.begin(), current.end(), [](char c) { return isdigit(c); })) {
                int n, k;
                std::istringstream iss(current);
                iss >> n >> k >> low_bound >> high_bound;
            }
            else {
                code.push_back(current);
            }
        }
    }

    return bencodes;
}

int main() {
    auto codes = code_list("testing/grassl.txt");
    
    for (auto &code : codes) {
        BMatrix code_mat;
        int low_bound, high_bound, n, k;
        std::tie(code_mat, low_bound, high_bound) = code;

        n = code_mat.m / 2;
        k = code_mat.m / 2 - code_mat.n;

        std::cout << n << " " << k << " " << low_bound << " " << high_bound << " - ";
        int dist = get_distance_with_parallelized_middle(code_mat, COMPUTE_TYPE { true, false, 16 });
        std::cout << dist << std::endl;

        if ((dist < low_bound || dist > high_bound) && (low_bound != 0 && high_bound != 0)) {
            std::cout << "Error: " << dist << " " << low_bound << " " << high_bound << std::endl;
            break;
        }
    }

    return 0;
}
