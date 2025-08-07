#include <iostream>
#include <algorithm>
#include <functional>
#include <vector>

#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "quantum_utilities.hpp"



int symplectic_weight(BVector v) {
    int weight = 0;
    for (int i = 0; i < v.no_buckets(); ++i)
        weight+= popcount( (v.vec[i] & x_mask) | ((v.vec[i] & z_mask) >> 1) );
    return weight;
}

bool is_css(BMatrix code) {
    bool is_css = true;
    
    for (int i = 0; i < code.n; ++i) {
        bool is_x = false;
        bool is_z = false;
        for (int j = 0; j < code.row(i).no_buckets(); ++j) {
            is_x|= code.row(i).vec[j] & x_mask;
            is_z|= code.row(i).vec[j] & z_mask;    
        }
        is_css&= !(is_x && is_z);
    }

    return is_css;
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


int non_stab_dist(BMatrix stab_mat, BMatrix closure_mat) {
    const int n = stab_mat.m;

    to_row_echelon(stab_mat);
    to_row_echelon(closure_mat);
    
    try {
        for (int d = 1; d <= n; ++d) {
            combinations(n, d, [&](std::vector<bool> vec) {
                BVector word(vec);
                if (!in_span(stab_mat, word) && in_span(closure_mat, word))
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


std::pair<int, int> bruteforce_zx_distance0(BMatrix stab_mat) {
    BMatrix closure_mat, z_stab, x_stab, z_closure, x_closure;

    closure_mat = isotropic_closure(stab_mat);
    std::tie(z_stab, x_stab) = zx_parts(stab_mat);
    std::tie(z_closure, x_closure) = zx_parts(closure_mat);

    return std::make_pair(
        non_stab_dist(z_stab, z_closure),
        non_stab_dist(x_stab, x_closure)
    );
}


int bruteforce_distance1(BMatrix stab_mat) {
    int z_dist, x_dist;

    std::tie(z_dist, x_dist) = bruteforce_zx_distance0(stab_mat);

    return std::min(z_dist, x_dist);
}


std::pair<BVector, BVector> zx_parts(BVector vec) {
    BVector z_vec(vec.n / 2), x_vec(vec.n / 2);
    for (int i = 0; i < vec.n; ++i) {
        z_vec.set(i, vec.get(2 * i));
        x_vec.set(i, vec.get(2 * i + 1));
    }

    return std::make_pair(z_vec, x_vec);
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


BMatrix logical_operators(BMatrix v_base) {
    BMatrix result, extension = basis_completion(v_base);

    while (!extension.empty()) {
        BVector last_ext = extension.last_row();

        int anticommuter = -1;
        for (int i = 0; i < v_base.n; ++i) {
            if (sym_prod(last_ext, v_base.row(i))) {
                anticommuter = i;
            }
        }

        if (anticommuter == -1) {
            result.append_row(last_ext); // if the operator commutes with everyone, we can just add it to the closure
        }
        else {
            // otherwise, if it does not commute with some anticommuter in v_base, we can use it to make every other vector in the extension commute with the anticommuter and then pop
            for (int i = 0; i < extension.n - 1; ++i) {
                if (sym_prod(extension.row(i), v_base.row(anticommuter)))
                    extension.add_rows(extension.n - 1, i);
            }
        }
        extension.pop_row();
    }

    for (int i = 0; i < result.n; i+= 2) {
        int anticommuter = -1;
        for (int j = i + 1; j < result.n; ++j) {
            if (sym_prod(result.row(i), result.row(j))) {
                anticommuter = j;
                break;
            }
        }
        if (anticommuter == -1) {
            std::cerr << "Logical operators are not anticommuting" << std::endl;
            my_assert(false);
        }

        std::swap(result.row(i + 1), result.row(anticommuter));
        for (int j = i + 2; j < result.n; ++j) {
            if (sym_prod(result.row(i), result.row(j)))
                result.add_rows(i + 1, j);
            if (sym_prod(result.row(i + 1), result.row(j)))
                result.add_rows(i, j);
        }
    }

    return result;
}

BMatrix destabilizers(BMatrix stab, BMatrix logs) {
    BMatrix comm_total;
    comm_total = stab;
    for (int i = 0; i < logs.n; ++i)
        comm_total.append_row(logs.row(i));

    BMatrix destab = basis_completion(comm_total);

    my_assert(destab.n == stab.n);

    // first, we ensure that the destabilizers commute with all logical operators
    for (int i = 0; i < destab.n; ++i) {
        for (int j = 0; j < logs.n; j+= 2) {
            if (sym_prod(destab.row(i), logs.row(j)))
                destab.row(i) = destab.row(i) + logs.row(j + 1);
            if (sym_prod(destab.row(i), logs.row(j + 1)))
                destab.row(i) = destab.row(i) + logs.row(j);
        }
    }

    // now we pair up a destabilizer to each stabilizer
    for (int i = 0; i < stab.n; ++i) {
        int anticomm = -1;
        for (int j = i; j < destab.n; ++j) {
            if (sym_prod(stab.row(i), destab.row(j))) {
                anticomm = j;
                break;
            }
        }
        if (anticomm == -1) {
            std::cerr << "Unpaired destabilizer" << std::endl;
            my_assert(false);
        }

        std::swap(destab.row(i), destab.row(anticomm));
        for (int j = i + 1; j < destab.n; ++j) {
            if (sym_prod(stab.row(i), destab.row(j)))
                destab.add_rows(i, j);
        }
    }

    // Enforce that the destabilizers commute with the stabilizers
    for (int i = 0; i < destab.n; ++i) {
        for (int j = i + 1; j < destab.n; ++j) {
            if (sym_prod(destab.row(i), destab.row(j))) {
                destab.row(i) = destab.row(i) + stab.row(j);
            }
        }
    }

    return destab;
}
