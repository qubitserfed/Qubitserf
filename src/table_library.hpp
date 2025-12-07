#pragma once

#include "linear_algebra.hpp"
#include <vector>
#include <utility>

/*
    CSS Code representation
    Important note:
        the z_stabs, x_stabs, z_logicals, x_logicals are not corresponding to the code provided as inputs to css_parse_code,
        but to an isomorphic code (column-permuted), with the pivots on the right. To recover a generator matrix for the 
        Z stabilizers for example, one must apply:
            apply_perm(z_stabs, invert_perm(z_perm));
*/

using Perm = std::vector<int>;

struct CSS_Code {
    BMatrix z_stabs, x_stabs;
    BMatrix z_logicals, x_logicals;
    std::vector<int> z_pivots, x_pivots;
    Perm z_perm, x_perm;
};

struct CSS_Table {
    std::vector<BVector> x_lookup, z_lookup;
};

using i64 = long long;

// Function Declarations
i64 comb(int n, int k);
Perm invert_perm(Perm p);
Perm pivots_right(int n, std::vector<int> pivots);
void apply_perm(BMatrix &mat, Perm perm);
void apply_perm(BVector &vec, Perm perm);
Perm pivots_to_right(int n, std::vector<int> pivots);
BVector parse_word(char *word, int n);

// Core functions
CSS_Code* css_parse_code(char **code, int n, int k);
CSS_Table* css_make_table(int n, int k, int bd, CSS_Code *code);
std::pair<int, int> css_lookup_word(BVector word, CSS_Code *code, CSS_Table *table);

