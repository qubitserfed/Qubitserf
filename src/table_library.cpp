#include "utility.hpp"
#include "quantum_utilities.hpp"
#include "linear_algebra.hpp"
#include "combinatorics.hpp"

#include <algorithm>
#include <cstdlib>
#include <set>

#include "table_library.hpp"



i64 comb(int n, int k) {
    i64 res = 1;
    for (int i = 1; i <= k; ++i)
        res = res * (n - i + 1) / i;
    return res;
}

Perm invert_perm(Perm p) {
    Perm res(p.size());
    for (int i = 0; i < (int)p.size(); ++i) {
        res[p[i]] = i;
    }
    return res;
}

/*
    Generates a permutation of 0..n-1 such that the elements of pivots are on the right, maintaining the relative order of the pivots
    (This is a stable partition function)
*/
Perm pivots_right(int n, std::vector<int> pivots) {
    Perm res(n);
    std::set<int> positions;
    for (int i = 0; i < n; ++i)
        positions.insert(i);

    for (int i = 0; i < (int)pivots.size(); ++i) {
        int pivot = pivots[int(pivots.size()) - i - 1];

        res[n - i - 1] = pivot;
        positions.erase(pivot);
    }
    
    for (int i = 0; i < n - (int)pivots.size(); ++i) {
        int pos = *positions.begin();
        positions.erase(pos);
        res[i] = pos;
    }

    return res;
}

void apply_perm(BMatrix &mat, Perm perm) {
    Perm inv_perm = invert_perm(perm);
    int num_rows = mat.n;
    int num_cols = mat.m;
    BMatrix temp_mat(num_rows, num_cols);

    for (int r = 0; r < num_rows; ++r) {
        for (int c = 0; c < num_cols; ++c) {
            temp_mat.set(r, c, mat.get(r, perm[c]));
        }
    }
    mat = temp_mat;
}

void apply_perm(BVector &vec, Perm perm) {
    int num_cols = vec.n;
    BVector temp_vec(num_cols);

    for (int c = 0; c < num_cols; ++c)
        temp_vec.set(c, vec.get(perm[c]));
    vec = temp_vec;
}

CSS_Code* css_parse_code(char **code, int n, int k) {
    BMatrix stab_mat;

    for (int i = 0; i < n - k; ++i) {
        BVector row(2 * n);

        for (int j = 0; j < n; ++j) {
            char c = code[i][j];
            if (c == 'I') {
                row.set(2 * j, false);
                row.set(2 * j + 1, false);
            }
            else if (c == 'X') {
                row.set(2 * j, false);
                row.set(2 * j + 1, true);
            }
            else if (c == 'Y') {
                row.set(2 * j, true);
                row.set(2 * j + 1, true);
            }
            else if (c == 'Z') {
                row.set(2 * j, true);
                row.set(2 * j + 1, false);
            }
            else {
                std::cerr << "Unrecognized character in stabilizer matrix: " << c << std::endl;
                throw std::invalid_argument("Unrecognized character in stabilizer matrix");
            }
        }

        stab_mat.append_row(row);
    }

    BMatrix logs = logical_operators(stab_mat);

    BMatrix z_stabs, x_stabs; std::tie(z_stabs, x_stabs) = zx_parts(stab_mat);
    BMatrix z_logs, x_logs; std::tie(z_logs, x_logs) = zx_parts(logs);

    to_row_echelon(z_stabs); z_stabs.remove_zeros();
    to_row_echelon(x_stabs); x_stabs.remove_zeros();

    std::vector<int> z_pivots, x_pivots;
    for (int i = 0; i < z_stabs.n; ++i) {
        int pivot = -1;
        for (int j = 0; j < z_stabs.m; ++j) {
            if (z_stabs.get(i, j)) {
                pivot = j;
                break;
            }
        }
        z_pivots.push_back(pivot);
    }
    for (int i = 0; i < x_stabs.n; ++i) {
        int pivot = -1;
        for (int j = 0; j < x_stabs.m; ++j) {
            if (x_stabs.get(i, j)) {
                pivot = j;
                break;
            }
        }
        x_pivots.push_back(pivot);
    }

    Perm z_perm = pivots_right(n, z_pivots);
    Perm x_perm = pivots_right(n, x_pivots);

    apply_perm(z_stabs, z_perm);
    apply_perm(x_stabs, x_perm);

    apply_perm(z_logs, z_perm);
    apply_perm(x_logs, x_perm);

    return new CSS_Code{
        z_stabs,
        x_stabs,
        z_logs,
        x_logs,
        z_pivots,
        x_pivots,
        z_perm,
        x_perm
    };
}

BVector parse_word(char *word, int n) {
    BVector vec(2 * n);
    for (int i = 0; i < n; ++i) {
        switch (word[i]) {
            case 'I':
                vec.set(2 * i, 0);
                vec.set(2 * i + 1, 0);
                break;
            case 'Z':
                vec.set(2 * i, 1);
                vec.set(2 * i + 1, 0);
                break;
            case 'X':
                vec.set(2 * i, 0);
                vec.set(2 * i + 1, 1);
                break;
            case 'Y':
                vec.set(2 * i, 1);
                vec.set(2 * i + 1, 1);
                break;
            default:
                std::cerr << "Unrecognized character in stabilizer string: " << word[i] << std::endl;
                throw std::invalid_argument("Unrecognized character in stabilizer string");
        }
    }
    return vec;
}

CSS_Table *css_make_table(int n, int k, int bd, CSS_Code *code) {
    CSS_Table *table = new CSS_Table;

    table->x_lookup.push_back(BVector(0));
    table->z_lookup.push_back(BVector(0));

    for (int d = 1; d <= bd; ++d) {
        table->x_lookup.push_back(BVector(1ULL << (n - code->z_stabs.n)));
        table->z_lookup.push_back(BVector(1ULL << (n - code->x_stabs.n)));

        parallel_combinations(n, d, [&](BVector &v) {
            BVector z_word = v, x_word = v;

            apply_perm(z_word, code->z_perm);
            apply_perm(x_word, code->x_perm);

            for (int i = 0; i < int(code->z_stabs.n); ++i)
                if (z_word.get(n - i - 1))
                    z_word = z_word + code->z_stabs.row(code->z_stabs.n - i - 1);
            z_word.resize(n - int(code->z_stabs.n));
            
            for (int i = 0; i < int(code->x_stabs.n); ++i)
                if (x_word.get(n - i - 1))
                    x_word = x_word + code->x_stabs.row(code->x_stabs.n - i - 1);
            x_word.resize(n - int(code->x_stabs.n));

            u64 lookup_index_z = z_word.n > 0 ? z_word.vec[0] : 0;
            u64 lookup_index_x = x_word.n > 0 ? x_word.vec[0] : 0;

            table->x_lookup.back().set(lookup_index_x, 1); //until we find what's broken
            table->z_lookup.back().set(lookup_index_z, 1);
            return false;

            // Thread-safe update using atomic fetch_or
            if (lookup_index_x < (u64)table->x_lookup.back().n) {
                u64 word_idx = lookup_index_x / 64;
                u64 bit_idx = lookup_index_x % 64;
                if (word_idx < table->x_lookup.back().vec.size()) {
                    __atomic_fetch_or(&table->x_lookup.back().vec[word_idx], (1ULL << bit_idx), __ATOMIC_RELAXED);
                }
            }
            
            if (lookup_index_z < (u64)table->z_lookup.back().n) {
                u64 word_idx = lookup_index_z / 64;
                u64 bit_idx = lookup_index_z % 64;
                if (word_idx < table->z_lookup.back().vec.size()) {
                    __atomic_fetch_or(&table->z_lookup.back().vec[word_idx], (1ULL << bit_idx), __ATOMIC_RELAXED);
                }
            }
            return false;
        }, 2048);
    }

    return table;
}

std::pair<int, int> css_lookup_word(BVector word, CSS_Code *code, CSS_Table *table) {
    const int n = code->z_stabs.m;
    BVector z_word, x_word; std::tie(z_word, x_word) = zx_parts(word);

    apply_perm(z_word, code->z_perm);
    apply_perm(x_word, code->x_perm);

    for (int i = 0; i < int(code->z_stabs.n); ++i)
        if (z_word.get(n - i - 1))
            z_word = z_word + code->z_stabs.row(code->z_stabs.n - i - 1);
    z_word.resize(n - int(code->z_stabs.n));
    
    for (int i = 0; i < int(code->x_stabs.n); ++i)
        if (x_word.get(n - i - 1))
            x_word = x_word + code->x_stabs.row(code->x_stabs.n - i - 1);
    x_word.resize(n - int(code->x_stabs.n));

    i64 lookup_index_x = x_word.n > 0 ? x_word.vec[0] : 0;
    i64 lookup_index_z = z_word.n > 0 ? z_word.vec[0] : 0;

    int ans_x = -1, ans_z = -1;
    
    if (z_word.is_zero())
        ans_z = 0;
    else {
        for (int d = 1; d < (int)table->z_lookup.size(); ++d) {
            if (table->z_lookup[d].get(lookup_index_z)) {
                ans_z = d;
                break;
            }
        }
    }

    if (x_word.is_zero())
        ans_x = 0;
    else {
        for (int d = 1; d < (int)table->x_lookup.size(); ++d) {
            if (table->x_lookup[d].get(lookup_index_x)) {
                ans_x = d;
                break;
            }
        }
    }

    return std::make_pair(ans_z, ans_x);
}

extern "C" void *css_get_code(char **code, int n, int k) {
    return css_parse_code(code, n, k);
}

extern "C" void *css_make_table(int n, int k, int bd, void *code) {
    return css_make_table(n, k, bd, code);
}

extern "C" void css_free_table(void *table) {
    delete table;
}

extern "C" void css_free_code(void *code) {
    delete code;
}

extern "C" int *css_lookup_word_str(char *stab, int n, void *code, void *table) {
    BVector vec = parse_word(stab, n);

    CSS_Code *code_ = (CSS_Code *)code;
    CSS_Table *table_ = (CSS_Table *)table;

    auto res = css_lookup_word(vec, code_, table_);

    return new int[] { res.first, res.second };
}

/*
Potential additional functionality:
    - storing lookup tables in memory

*/