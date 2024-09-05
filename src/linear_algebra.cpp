#include <algorithm>
#include <iostream>
#include <vector>

#include <cassert>

#include "linear_algebra.hpp"

using u64 = unsigned long long;

void my_assert(bool arg) {
    assert(arg);
}

int popcount(u64 num) {
    const u64 masks[] = {
        0x5555555555555555ULL,
        0x3333333333333333ULL,
        0xf0f0f0f0f0f0f0fULL,
        0xff00ff00ff00ffULL,
        0xffff0000ffffULL,
        0xffffffffULL,
    };
    const int n = sizeof(masks) / sizeof(masks[0]);

    num = (num & masks[0]) + ((num >> (1 << 0)) & masks[0]);
    num = (num & masks[1]) + ((num >> (1 << 1)) & masks[1]);
    num = (num & masks[2]) + ((num >> (1 << 2)) & masks[2]);
    num = (num & masks[3]) + ((num >> (1 << 3)) & masks[3]);
    num = (num & masks[4]) + ((num >> (1 << 4)) & masks[4]);
    num = (num & masks[5]) + ((num >> (1 << 5)) & masks[5]);

    return num;
}

/// Start of BVector definitions

int BVector::no_buckets() {
    return vec.size();
}

bool BVector::get(int pos) {
    return vec[pos / 64] & (1ULL << (pos % 64));
}

bool BVector::is_zero() {
    bool anyone_nonzero = false;
    for (const auto &ll: vec)
        anyone_nonzero|= ll;
    return !anyone_nonzero;
}

void BVector::set(int pos, bool val) {
    bool old_val = get(pos);
    if (old_val != val)
        vec[pos / 64]^= (1ULL << (pos % 64));
}

int BVector::weight() {
    int result = 0;
    for (u64 entry: vec)
        result+= popcount(entry);
    return result;
}

BVector::BVector() : n(0) { }

BVector::BVector(int _n) : n(_n) {
    vec.resize((n + 63) / 64);
}

BVector::BVector(std::vector<bool> v) {
    n = v.size();
    vec.resize((n + 63) / 64);
    for (int i = 0; i < n; ++i)
        set(i, v[i]);
}

void print(BVector vec) {
    const int n = vec.n;

    for (int i = 0; i < n; ++i)
        std::cout << vec.get(i) << " \n"[i == n - 1];
}

BVector operator + (const BVector &a, const BVector &b) {
    my_assert(a.n == b.n);

    const int n = a.n;
    BVector res(n);
    for (int i = 0; i < a.vec.size(); ++i)
        res.vec[i] = (a.vec[i] ^ b.vec[i]);
    return res;
}

bool operator * (const BVector &a, const BVector &b) {
    my_assert(a.n == b.n);

    const int n = a.n;
    int ans = 0;
    for (int i = 0; i < a.vec.size(); ++i)
        ans+= popcount(a.vec[i] & b.vec[i]);
    
    return ans % 2;
}

/// End of BVector definitions
/// -------------------------------------------------------
/// Start of BMatrix definitions


BVector &BMatrix::row(int i) {
    my_assert(0 <= i && i < n);

    return mat[i];
}

BVector &BMatrix::last_row() {
    my_assert(n > 0);
    
    return mat.back();
}

void BMatrix::append_row(BVector v) {
    my_assert(v.n == m || mat.empty());

    if (mat.empty())
        m = v.n;

    mat.push_back(v);
    n+= 1;
}

void BMatrix::pop_row() {
    mat.pop_back();
    n-= 1;
}

void BMatrix::set(int i, int j, bool val) {
    mat[i].set(j, val);
}

bool BMatrix::get(int i, int j) {
    return mat[i].get(j);
}

bool BMatrix::empty() {
    return mat.empty();
}

void BMatrix::add_rows(int i, int j) { // add row i to row j
    for (int k = 0; k < mat[i].no_buckets(); ++k) {
        mat[j].vec[k]^= mat[i].vec[k];
    }
}

void BMatrix::swap_rows(int i, int j) {
    std::swap(mat[i], mat[j]);
}

void BMatrix::swap_cols(int i, int j) {
    for (int k = 0; k < mat.size(); ++k) {
        const bool a = get(k, i);
        const bool b = get(k, j);
        set(k, i, b);
        set(k, i, a);
    }
}

void BMatrix::remove_zeros() {
    for (int i = 0; i < int(mat.size()); ++i) {
        if (mat[i].weight() == 0) {
            std::swap(mat[i], mat.back());
            mat.pop_back();
            i-= 1;
        }
    }

    n = mat.size();
    if (n == 0)
        m = 0;
}

void BMatrix::sort_rows() {
    std::sort(mat.begin(), mat.end(), [&](const BVector &a, const BVector &b){
        return a.vec > b.vec;
    });
}

BMatrix::BMatrix() : n(0), m(0) { }

BMatrix::BMatrix(int _n, int _m) : n(_n), m(_m) {
    mat.resize(n);
    for (auto &row: mat)
        row = BVector(m);
}


void print(BMatrix mat) {
    const int n = mat.n;
    const int m = mat.m;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j)
            std::cout << mat.get(i, j) << " \n"[j == m - 1];
    }
    std::cout << std::endl;
}

BMatrix identity(const int &n) {
    BMatrix res(n, n);
    for (int i = 0; i < n; ++i)
        res.set(i, i, 1);
    return res;
}

BMatrix transpose(BMatrix &origin) {
    const int n_origin = origin.n;
    const int m_origin = origin.m;
    const int n = m_origin;
    const int m = n_origin;

    BMatrix res(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            res.set(i, j, origin.get(j, i));

    return res;
}

BVector operator * (BVector vec, BMatrix mat) {
    my_assert(vec.n == mat.n);

    const int n = mat.n;
    const int m = mat.m;
    
    BVector res(m);
    for (int j = 0; j < m; ++j) {
        bool rj = 0;
        for (int i = 0; i < n; ++i)
            rj^= (vec.get(i) & mat.get(i, j));
        res.set(j, rj);
    } 
    return res;
}

/// End of Bmatrix definitions
/// -------------------------------------------------------
/// Start of Linear algebra utilities

// input: a matrix
// output: the rank of the matrix
// effect: the matrix gets put into row_echelon form solely via elementary row operations
// invariant: if the input is CSS, the output will be CSS as well
int to_row_echelon(BMatrix &mat) {
    const int n = mat.n;
    const int m = mat.m;

    int rank = 0;
    int pivot_row = -1;
    int pivot_col = -1; // at the start of each loop, the submatrix [0..pivot_row][0..pivot_col] in in row_echelon form
    while(pivot_row < n - 1 && pivot_col < m - 1) {
        bool found_nonzero = false;
        
        for (int i = pivot_row + 1; i < n; ++i) {
            if (mat.get(i, pivot_col + 1) == 1) {
                found_nonzero = true;
                mat.swap_rows(pivot_row + 1, i);
                pivot_row+= 1;
                pivot_col+= 1;
                break;
            }
        }

        if (found_nonzero) {
            rank+= 1;
            for (int i = 0; i < n; ++i) {
                if (i != pivot_row && mat.get(i, pivot_col))
                    mat.add_rows(pivot_row, i);
            }
        }
        else {
            pivot_col+= 1;
        }
    }

    return rank;
}

// input:  a matrix M and a list of columns L
// effect: the submatrix given by the columns in L is put into row-echelon form
// output: a sublist of L consisting of the columns where the first nonzero entry of the rows appear
std::vector<int> restricted_row_echelon(BMatrix &mat, std::vector<int> cols) {
    const int n = mat.n;
    const int m = mat.m;

    int rank = 0;
    int pivot_row = -1;
    int pivot_col = -1;
    std::vector<int> returned_columns;

    while(pivot_row < n - 1 && pivot_col < int(cols.size()) - 1) {
        bool found_nonzero = false;
        
        for (int i = pivot_row + 1; i < n; ++i) {
            if (mat.get(i, cols[pivot_col + 1]) == 1) {
                found_nonzero = true;
                mat.swap_rows(pivot_row + 1, i);
                pivot_row+= 1;
                pivot_col+= 1;
                break;
            }
        }

        if (found_nonzero) {
            returned_columns.push_back(cols[pivot_col]);
            for (int i = 0; i < n; ++i) {
                if (i != pivot_row && mat.get(i, cols[pivot_col]))
                    mat.add_rows(pivot_row, i);
            }
        }
        else {
            pivot_col+= 1;
        }
    }

    return returned_columns;
}


// input: an n-dimensional vector v and an nxm matrix M in row_echelon form
// output: a canonical representative of the coset v + [the space spanned by the the rows of M]
// explainer: by canonical I mean that if V = [the space spanned by the the rows of M] and for two
// vectors v and w, v + V = w + V, canonical_quotient(v, mat) will equal canonical_quotient(w, mat) 
BVector canonical_quotient(BVector vec, BMatrix &mat) {
    my_assert(vec.n == mat.m);

    const int n = mat.n;
    const int m = mat.m;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (mat.get(i, j)) {
                if (vec.get(j))
                    vec = vec + mat.row(i);
                break;
            }
        }
    }

    return vec;
}

// input: an n-dimensional vector v and an nxm matrix M in row_echelon form
// output: whether v is part of the subspace spanned by the rows of M
bool in_span(BMatrix &mat, const BVector &v) {
    return canonical_quotient(v, mat).is_zero();
}

// input: two matrices A and B
// output: the product A * B^T
BMatrix transposed_product(BMatrix &a, BMatrix &b) {
    my_assert(a.m == b.m);

    const int n = a.n;
    const int m = b.m;

    BMatrix res(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            res.set(i, j, a.row(i) * b.row(j));
    
    return res;
}

// input: a vector v and a matrix M
// output: the product v * M^T
BVector transposed_product(const BVector &v, BMatrix &mat) {
    my_assert(v.n == mat.m);

    const int n = mat.n;
    const int m = mat.m;

    BVector res(n);
    for (int i = 0; i < n; ++i)
        res.set(i, v * mat.row(i));

    return res;
}

// input: a matrix whose rows span a space V
// output: a matrix whose rows span a space W such that V and W
// intersect trivially and their sum is the total space
// invariant: if the input is CSS, the output will be CSS as well
BMatrix basis_completion(BMatrix mat) {
    const int n = mat.n;
    const int m = mat.m;

    to_row_echelon(mat);

    BMatrix res;
    for (int i = 0; i < m; ++i) {
        BVector candidate(m);

        candidate.set(i, true);
        candidate = canonical_quotient(candidate, mat);
        if (!candidate.is_zero()) {
            res.append_row(candidate);

            mat.append_row(candidate); // this needs to be optimized eventually
            to_row_echelon(mat); 
        }
    }

    return res;
}

constexpr u64 raw_sym_prod_ll(const u64 &a, const u64 &b) {
//  xz
    const u64 prods = ((a & z_mask) & ((b & x_mask) >> 1)) ^ ((b & z_mask) & ((a & x_mask) >> 1));
    return prods;
}

// input: two vectors of the same *even length*
// output: the value of the symplectic form over a and b (a \lambda_n b^T)
bool sym_prod(BVector a, BVector b) {
    my_assert(a.n == b.n);
    my_assert(a.n % 2 == 0);

    u64 acc = 0;
    for (int i = 0; i < a.no_buckets(); ++i)
        acc^= raw_sym_prod_ll(a.vec[i], b.vec[i]);
    return popcount(acc) % 2;
}

// input:  a matrix in whose rows span an isotropic space V
// output: a matrix whose rows span the isotropic closure of V
// warning: slowness inherited from basis_completion
BMatrix isotropic_closure(BMatrix res) {
    to_row_echelon(res);
    BMatrix extension = basis_completion(res);

    while (!extension.empty()) {
        BVector last_ext = extension.last_row();
        
        int anticommuter = -1;
        for (int i = 0; i < res.n; ++i) {
            if (sym_prod(last_ext, res.row(i))) {
                anticommuter = i;
            }
        }

        if (anticommuter == -1) {
            // if the operator commutes with everyone, we can just add it to the closure
            res.append_row(last_ext);
            extension.pop_row();
        }
        else {
            // otherwise, if it anticommutes with some anticommuter, we can use it to make every other vector in the extension commute with the anticommuter and then pop
            for (int i = 0; i < extension.n - 1; ++i) {
                if (sym_prod(extension.row(i), res.row(anticommuter)))
                    extension.add_rows(extension.n - 1, i);
            }
            extension.pop_row();
        }
    }

    to_row_echelon(res);

    return res;
}

// input: a matrix whose rows span a space V and a variable 'prioritary'
// output: a matrix whose rows span a space W such that V and W
// intersect trivially and their sum is the total space. If
// prioritary equals false, it will prioritize adding even-column rows,
// otherwise it will prioritize odd-column rows
// invariant: if the input is CSS, the output will be CSS as well
BMatrix prefered_basis_completion(BMatrix mat, bool prioritary) {
    const int n = mat.n;
    const int m = mat.m;

    to_row_echelon(mat);

    BMatrix res;
    for (int i = int(prioritary); i < m; i+= 2) {
        BVector candidate(m);

        candidate.set(i, true);
        candidate = canonical_quotient(candidate, mat);
        if (!candidate.is_zero()) {
            res.append_row(candidate);

            mat.append_row(candidate); // this needs to be optimized eventually
            to_row_echelon(mat); 
        }
    }

    for (int i = int(!prioritary); i < m; i+= 2) {
        BVector candidate(m);

        candidate.set(i, true);
        candidate = canonical_quotient(candidate, mat);
        if (!candidate.is_zero()) {
            res.append_row(candidate);

            mat.append_row(candidate); // this needs to be optimized eventually
            to_row_echelon(mat); 
        }
    }

    return res;
}


BMatrix prefered_isotropic_closure(BMatrix res, bool prioritary) {
    to_row_echelon(res);
    BMatrix extension = prefered_basis_completion(res, prioritary);

    std::reverse(extension.mat.begin(), extension.mat.end());

    while (!extension.empty()) {
        BVector last_ext = extension.last_row();
        
        int anticommuter = -1;
        for (int i = 0; i < res.n; ++i) {
            if (sym_prod(last_ext, res.row(i))) {
                anticommuter = i;
            }
        }

        if (anticommuter == -1) {
            // if the operator commutes with everyone, we can just add it to the closure
            res.append_row(last_ext);
            //print(res);
            extension.pop_row();
        }
        else {
            // otherwise, if it anticommutes with some anticommuter, we can use it to make every other vector in the extension commute with the anticommuter and then pop
            for (int i = 0; i < extension.n - 1; ++i) {
                if (sym_prod(extension.row(i), res.row(anticommuter)))
                    extension.add_rows(extension.n - 1, i);
            }
            extension.pop_row();
        }
    }

    to_row_echelon(res);

    return res;
}
