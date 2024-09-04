#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <cassert>

using u64 = unsigned long long;

const u64 x_mask = 0xAAAAAAAAAAAAAAAAULL;
const u64 z_mask = ~x_mask;

void my_assert(bool arg) {
    assert(arg);
}

int popcount(u64 num) { // unoptimized piece of crap (junk, fix this)
    int ans = 0;
    for (int i = 0; i < 64; ++i)
        ans+= (num & (1ULL << i)) ? 1 : 0;
    return ans;
}

/// Start of BVector definitions

struct BVector {
    std::vector<u64> vec;
    int n;

    int no_buckets() {
        return vec.size();
    }

    bool get(int pos) {
        return vec[pos / 64] & (1ULL << (pos % 64));
    }

    bool is_zero() {
        bool anyone_nonzero = false;
        for (const auto &ll: vec)
            anyone_nonzero|= ll;
        return !anyone_nonzero;
    }

    void set(int pos, bool val) {
        bool old_val = get(pos);
        if (old_val != val)
            vec[pos / 64]^= (1ULL << (pos % 64));
    }

    int weight() {
        int result = 0;
        for (u64 entry: vec)
            result+= popcount(entry);
        return result;
    }

    BVector() : n(0) { }

    BVector(int _n) : n(_n) {
        vec.resize((n + 63) / 64);
    }

    BVector(std::vector<bool> v) {
        n = v.size();
        vec.resize((n + 63) / 64);
        for (int i = 0; i < n; ++i)
            set(i, v[i]);
    }
};

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

struct BMatrix {
    std::vector<BVector> mat;
    int n, m;

    BVector &row(int i) {
        my_assert(0 <= i && i < n);

        return mat[i];
    }

    BVector &last_row() {
        my_assert(n > 0);
        
        return mat.back();
    }

    void append_row(BVector v) {
        my_assert(v.n == m || mat.empty());

        if (mat.empty())
            m = v.n;

        mat.push_back(v);
        n+= 1;
    }

    void pop_row() {
        mat.pop_back();
        n-= 1;
    }

    void set(int i, int j, bool val) {
        mat[i].set(j, val);
    }

    bool get(int i, int j) {
        return mat[i].get(j);
    }

    bool empty() {
        return mat.empty();
    }

    void add_rows(int i, int j) { // add row i to row j
        for (int k = 0; k < mat[i].no_buckets(); ++k) {
            mat[j].vec[k]^= mat[i].vec[k];
        }
    }

    void swap_rows(int i, int j) {
        std::swap(mat[i], mat[j]);
    }

    void swap_cols(int i, int j) {
        for (int k = 0; k < mat.size(); ++k) {
            const bool a = get(k, i);
            const bool b = get(k, j);
            set(k, i, b);
            set(k, i, a);
        }
    }

    void remove_zeros() {
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

    void sort_rows() {
        std::sort(mat.begin(), mat.end(), [&](const BVector &a, const BVector &b){
            return a.vec > b.vec;
        });
    }

    BMatrix() : n(0), m(0) { }
    BMatrix(int _n, int _m) : n(_n), m(_m) {
        mat.resize(n);
        for (auto &row: mat)
            row = BVector(m);
    }
};

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

BMatrix transpose(const BMatrix &origin) {
    const int n_origin = origin.n;
    const int m_origin = origin.m;
    const int n = m_origin;
    const int m = n_origin;

    BMatrix res(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            res.set(i, j, res.get(j, i));

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


/// End of Linear algebra utilities
/// -------------------------------------------------------
/// Start of Combinatorial utilities

// input:  a number n and a function whose sole argument is a vector of bools
// effect: the function gets called on all length n bitstrings in lexicographical order
void partitions(int n, std::function<void(std::vector<bool>)> f) {
    std::vector<bool> stack;
    stack.reserve(n);

    std::function<void()> bkt = [&]() {
        if (stack.size() == n) {
            f(stack);
            return;
        }

        stack.push_back(0);
        bkt();
        stack.pop_back();

        stack.push_back(1);
        bkt();
        stack.pop_back();
    };

    bkt();
}

// input:  a number n and a function whose sole argument is a vector of bools
// effect: the function gets called on all length n bitstrings of Hamming weight k in lexicographical order
void combinations(int n, int k, std::function<void(std::vector<bool>) > f) {
    std::vector<bool> stack;
    stack.reserve(n);

    std::function<void(int)> bkt = [&](int weight) {
        if (stack.size() == n) {
            f(stack);
            return;
        }

        if (weight <= n - (int(stack.size()) + 1)) {
            stack.push_back(0);
            bkt(weight);
            stack.pop_back();
        }

        if (weight > 0) {
            stack.push_back(1);
            bkt(weight - 1);
            stack.pop_back();
        }
    };

    bkt(k);
}

/// End of Combinatorial utilities
/// -------------------------------------------------------
/// Start of Quantum Error Correction Utilities

int stabilizer_weight(BVector v) {
    int weight = 0;
    for (int i = 0; i < v.no_buckets(); ++i)
        weight+= popcount( (v.vec[i] & x_mask) | ((v.vec[i] & z_mask) >> 1) );
    return weight;
}


// input: a number n and a function whose sole argument is a vector of bools
// effect: the function gets called on all length 2n bitstrings corresponding to weight k stabilizers in lexicographical order
void iterate_words_of_weight(int n, int k, std::function<void(std::vector<bool>) > f) {
    std::vector<bool> stack;
    stack.reserve(n);

    std::function<void(int)> bkt = [&](int weight) {
        if (stack.size() == 2 * n) {
            f(stack);
            return;
        }

        if (weight <= n - (int(stack.size()) / 2 + 1)) {
            stack.push_back(0); stack.push_back(0);
            bkt(weight);
            stack.pop_back(); stack.pop_back();
        }

        if (weight > 0) {
            stack.push_back(0); stack.push_back(1);
            bkt(weight - 1);
            stack.pop_back(); stack.pop_back();

            stack.push_back(1); stack.push_back(0);
            bkt(weight - 1);
            stack.pop_back(); stack.pop_back();

            stack.push_back(1); stack.push_back(1);
            bkt(weight - 1);
            stack.pop_back(); stack.pop_back();
        }
    };

    bkt(k);
}

int bruteforce_distance0(BMatrix gen_mat) {
    const int n = gen_mat.n;
    const int m = gen_mat.m;
    const int no_qubits = m / 2;

    to_row_echelon(gen_mat);

    BMatrix logical_ops = isotropic_closure(gen_mat);
    to_row_echelon(logical_ops);
    
    try {
        for (int d = 1; d <= no_qubits; ++d) {
            iterate_words_of_weight(no_qubits, d, [&](std::vector<bool> vec) {
                BVector word(vec);
                if (!in_span(gen_mat, word) && in_span(logical_ops, word))
                    throw d;
            });
        }
    }
    catch (int d) {
        return d;
    }

    my_assert(false);
    return 0;
}

std::pair<BMatrix, BMatrix> zx_parts(BMatrix mat) {
    const int n = mat.n;
    const int m = mat.m;
    const int no_qubits = mat.m / 2;

    BMatrix z_mat, x_mat;

    for (int i = 0; i < n; ++i) {
        std::vector<bool> z_vec, x_vec;
        
        for (int j = 0; j < no_qubits; ++j) {
            z_vec.push_back(mat.get(i, 2 * j));
            x_vec.push_back(mat.get(i, 2 * j + 1));
        }

        z_mat.append_row(z_vec);
        x_mat.append_row(x_vec);
    }

    return std::make_pair(z_mat, x_mat);
}

/// End of Quantum Error Correction Utilities
/// -------------------------------------------------------
/// Start of Brouwer-Zimmerman Algorithm

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

    for (int d = 1; d <= n; ++d) {
        for (auto gen_pair: gamma_seq) {
            BMatrix gamma_mat = gen_pair.first;

            combinations(n, d, [&](std::vector<bool> v0) {
                BVector vec(v0);

                BVector codeword = vec * gamma_mat;

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

int get_distance(BMatrix stab_mat) {
    BMatrix closed_mat, x_stab, z_stab, x_closed, z_closed;

    to_row_echelon(stab_mat);
    closed_mat = isotropic_closure(stab_mat);
    to_row_echelon(stab_mat);
    
    std::tie(z_stab, x_stab) = zx_parts(stab_mat);
    std::tie(z_closed, x_closed) = zx_parts(closed_mat);

    z_stab.remove_zeros(); z_stab.sort_rows();
    z_closed.remove_zeros(); z_closed.sort_rows();

    x_stab.remove_zeros(); x_stab.sort_rows();
    x_closed.remove_zeros(); x_closed.sort_rows();

//    print(z_stab);
//    print(z_closed);

    const int x_dist = brouwer_zimmerman(x_stab, x_closed);
    const int z_dist = brouwer_zimmerman(z_stab, z_closed);

    return std::min(z_dist, x_dist);
}


/// End of Brouwer-Zimmerman Algorithm
/// -------------------------------------------------------
/// Start of Testing

BMatrix steane_code() {
    const int n = 6;
    const int m = 7;

    // 0 is I, 1 is Z, 2 is X, 3 is Y
    const int entries[n][m] = {
        {0, 0, 0, 1, 1, 1, 1},
        {0, 0, 0, 3, 3, 3, 3},
        {0, 1, 1, 0, 0, 1, 1},
        {0, 3, 3, 0, 0, 3, 3},
        {1, 0, 1, 0, 1, 0, 1},
        {3, 0, 3, 0, 3, 0, 3}
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

BMatrix bmatrix_ben_conversion(std::vector<std::string> code) {
    const int m = 2 * code.back().size();
    const int n = code.size();

    BMatrix result;

    for (int i = 0; i < n; ++i) {
        BVector row(m);

        for (int j = 0; j < code[i].size(); ++j) {
            switch (code[i][j]) {
                case '_': break;
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

BMatrix bmatrix_standard_conversion(std::vector<std::string> code) {
    const int m = 2 * code.back().size();
    const int n = code.size();

    BMatrix result;

    for (int i = 0; i < n; ++i) {
        BVector row(m);

        for (int j = 0; j < code[i].size(); ++j) {
            switch (code[i][j]) {
                case 'I': break;
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

std::vector<BMatrix> ben_codes() {
    std::ifstream fi("d4_codes.txt");

    std::vector<BMatrix> bencodes;
    std::vector<std::string> code;
    std::string current;

    while (getline(fi, current)) {
        if (current == "") {
            bencodes.push_back(bmatrix_ben_conversion(code));
            code.clear();
        }
        else {
            code.push_back(current);
        }
    }

    return bencodes;
}


int main() {
    std::string current;
    std::vector<std::string> code_strs;
    
    while (getline(std::cin, current)) {
        if (current == "")
            continue;
        code_strs.push_back(current);
    }

    BMatrix code = bmatrix_standard_conversion(code_strs);

    std::cout << get_distance(code) << std::endl;

    return 0;
}
