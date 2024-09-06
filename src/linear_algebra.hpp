#pragma once

#include <algorithm>
#include <vector>
#include <cassert>

using u64 = unsigned long long;

const u64 x_mask = 0xAAAAAAAAAAAAAAAAULL;
const u64 z_mask = ~x_mask;

struct BVector {
    std::vector<u64> vec;
    int n;

    int no_buckets();
    bool get(int);
    bool is_zero();
    void set(int , bool);
    int weight();

    BVector();
    BVector(int);
    BVector(std::vector<bool>);
};

struct BMatrix {
    std::vector<BVector> mat;
    int n, m;

    BVector &row(int);

    BVector &last_row();

    void append_row(BVector);
    void pop_row();
    void set(int, int, bool);
    bool get(int, int);
    bool empty();

    void add_rows(int, int);
    void swap_rows(int, int);
    void swap_cols(int, int);
    void remove_zeros();
    void sort_rows();
    BMatrix();
    BMatrix(int, int);
};


void                my_assert                       (bool arg);
int                 popcount                        (u64 num);
void                print                           (BVector vec);
BVector             operator +                      (const BVector &, const BVector &);
bool                operator *                      (const BVector &, const BVector &);
void                print                           (BMatrix);
BMatrix             identity                        (const int &);
BMatrix             transpose                       (BMatrix &);
BVector             operator *                      (BVector, BMatrix);
int                 to_row_echelon                  (BMatrix &);
std::vector<int>    restricted_row_echelon          (BMatrix &, std::vector<int>);
BVector             canonical_quotient              (BVector, BMatrix &);
bool                in_span                         (BMatrix &, const BVector &);
BMatrix             transposed_product              (BMatrix &, BMatrix &);
BVector             transposed_product              (const BVector &, BMatrix &);
BMatrix             basis_completion                (BMatrix);
constexpr u64       raw_sym_prod_ll                 (const u64 &, const u64&);
bool                sym_prod                        (BVector, BVector);
BMatrix             isotropic_closure               (BMatrix);
