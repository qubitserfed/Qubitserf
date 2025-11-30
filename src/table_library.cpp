#include "utility.hpp"
#include "quantum_utilities.hpp"
#include "linear_algebra.hpp"

struct CSS_table {
    std::vector<ParallelHashTable*> ht_X;
    std::vector<ParallelHashTable*> ht_Z;
};

BMatrix parse_code(char **code, int n, int k) {
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

    return stab_mat;
}


extern "C" void *create_hashtable_css(int n, int k, int bd, char **stabs) {
    BMatrix code = parse_code(stabs, n, k);
    BMatrix logicals = logical_operators(code);

    BMatrix z_logicals, x_logicals; std::tie(z_logicals, x_logicals) = zx_parts(logicals);
    BMatrix z_stabilizers, x_stabilizers; std::tie(z_stabilizers, x_stabilizers) = zx_parts(code);

    std::vector<ParallelHashTable*> xhashes, zhashes;

    for (int weight = 0; weight <= bd; ++weight) {
        ParallelHashTable* ht_X = new ParallelHashTable(); ht_X->reset_for_code(n, bd, true);
        ParallelHashTable* ht_Z = new ParallelHashTable(); ht_Z->reset_for_code(n, bd, true);
        
        parallel_combinations(n, bd, [&](BVector &comb) {
            BVector anticomm_part, log_part;
            
            // X-part
            anticomm_part = comb * z_stabilizers;
            log_part = comb * z_logicals;
            ht_X->insert(anticomm_part, log_part);

            // Z-part
            anticomm_part = comb * x_stabilizers;
            log_part = comb * x_logicals;
            ht_Z->insert(anticomm_part, log_part);
            return false;
        }, 1024);
        
        xhashes.push_back(ht_X);
        zhashes.push_back(ht_Z);
    }

    return new CSS_table{xhashes, zhashes};
}

extern "C" int *query_hashtable_css(void *table_ptr, int n, int k, char *op) {
    CSS_table* table = static_cast<CSS_table*>(table_ptr);
    BVector op_vec(2 * n);

    for (int j = 0; j < n; ++j) {
        char c = op[j];
        if (c == 'I') {
            op_vec.set(2 * j, false);
            op_vec.set(2 * j + 1, false);
        }
        else if (c == 'X') {
            op_vec.set(2 * j, false);
            op_vec.set(2 * j + 1, true);
        }
        else if (c == 'Y') {
            op_vec.set(2 * j, true);
            op_vec.set(2 * j + 1, true);
        }
        else if (c == 'Z') {
            op_vec.set(2 * j, true);
            op_vec.set(2 * j + 1, false);
        }
        else {
            std::cerr << "Unrecognized character in operator: " << c << std::endl;
            throw std::invalid_argument("Unrecognized character in operator");
        }
    }

    BVector anticomm_part, log_part;
    std::tie(anticomm_part, log_part) = zx_parts(op_vec);

    int dist_z = -1, dist_x = -1;
    for (int weight = 0; weight < table->ht_Z.size(); ++weight) {
        if (table->ht_Z[weight - 1]->find(anticomm_part, log_part)) {
            dist_z = weight;
            break;
        }
    }

    for (int weight = 0; weight < table->ht_X.size(); ++weight) {
        if (table->ht_X[weight - 1]->find(anticomm_part, log_part)) {
            dist_x = weight;
            break;
        }
    }

    return (int[]){dist_z, dist_x};
}