#include <algorithm>
#include <functional>
#include <vector>

#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "quantum_utilities.hpp"



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

        if (weight <= (2 * n - int(stack.size())) / 2 - 1 ) {
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
