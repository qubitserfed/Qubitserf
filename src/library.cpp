#include <iostream>
#include <tuple>
#include <utility>

#include "linear_algebra.hpp"
#include "brouwer_zimmerman.hpp"

extern "C" {
    // input: n, m (the matrix dimensions), compute_type (0 singlethread, 1 multithread, 2 GPU (not yet implemented))
    // output: output to a pointer to the start of an int array containing the z distance and the x distance
    int *get_zx_distances(int, int, int, char **);

    // input: n, m (the matrix dimensions), compute_type (0 singlethread, 1 multithread, 2 GPU (not yet implemented))
    // output: the distance of the code
    int get_distance(int, int, int, char **);


    int *get_zx_distances_raw(
        int,    // number of qubits
        int,    // number of z stabilizers
        u64 **, // z stabilizers
        int,    // number of x stabilizers
        u64 **, // x stabilizers
        int     // compute type (0 singlethread, 1 multithread, 2 GPU (not yet implemented))
    );

    int get_distance_raw(
        int,    // number of qubits
        int,    // number of z stabilizers
        u64 **, // z stabilizers
        int,    // number of x stabilizers
        u64 **, // x stabilizers
        int     // compute type (0 singlethread, 1 multithread, 2 GPU (not yet implemented))
    );
}

BMatrix from_c_array(int n, int m, char **stab_mat) {
    BMatrix mat(n, 2 * m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            switch (stab_mat[i][j]) {
                case 'I':
                    break;
                case 'Z':
                    mat.set(i, 2 * j, 1);
                    break;
                case 'X':
                    mat.set(i, 2 * j + 1, 1);
                    break;
                default:
                    my_assert(false);
            }
        }
    }
    return mat;
}

BMatrix from_raw(int no_qubits, int nz, u64 **z_stab, int nx, u64 **x_stab) {
    BMatrix mat;

    std::cout << "Raw Z-stabs:\n";
    for (int i = 0; i < nz; ++i) {
        for (int j = 0; j < (no_qubits + 63) / 64; ++j)
            std :: cout << z_stab[i][j] << " ";
        std :: cout << '\n';
    }
    std::cout << "Raw X-stabs:\n";
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < (no_qubits + 63) / 64; ++j)
            std :: cout << x_stab[i][j] << " ";
        std :: cout << '\n';
    }

    std::cout << "Z stabilizers:\n";
    for (int i = 0; i < nz; ++i) {
        BVector new_stab(2*no_qubits);
        for (int j = 0; j < no_qubits; ++j)  {
            const int bucket = j / 64;
            const int bit = j % 64;

            int j2 = j * 2;
            const int bucket2 = j2 / 64;
            const int bit2 = j2 % 64;

            std :: cout << z_stab[i][bucket] << ' ' << j << ' ' << bit << ' ' << bucket << ' ' << (z_stab[i][bucket] & (1ULL << bit)) << std::endl;

            if (z_stab[i][bucket] & (1ULL << bit))
                new_stab.vec[bucket2]|= (1ULL << bit2);
        }
        mat.append_row(new_stab);
    }


    std::cout << "X stabilizers:\n";
    for (int i = 0; i < nx; ++i) {
        BVector new_stab(2*no_qubits);
        for (int j = 0; j < no_qubits; ++j)  {
            const int bucket = j / 64;
            const int bit = j % 64;

            new_stab.set(2*j+1, x_stab[i][bucket] & (1ULL << bit));
        }
        mat.append_row(new_stab);
    }

    std::cout << "Whole matrix:\n";

    return mat;
}

int *get_zx_distances_raw(int no_qubits, int no_z_stab, u64 **z_stabs_raw, int no_x_stab, u64 **x_stabs_raw, int compute_type) {
    BMatrix stab_mat = from_raw(no_qubits, no_z_stab, z_stabs_raw, no_x_stab, x_stabs_raw);
    int *res = new int[2];
    std::tie(res[0], res[1]) = get_zx_distances(stab_mat, (COMPUTE_TYPE)compute_type);
    return res;
}

int get_distance_raw(int no_qubits, int no_z_stab, u64 **z_stabs_raw, int no_x_stab, u64 **x_stabs_raw, int compute_type) {
    return get_distance(
        from_raw(no_qubits, no_z_stab, z_stabs_raw, no_x_stab, x_stabs_raw),
        (COMPUTE_TYPE) compute_type
    );
}

int *get_zx_distances(int n, int m, int compute_type, char **stab_mat) {
    BMatrix mat = from_c_array(n, m, stab_mat);
    int *res = new int[2];
    std::tie(res[0], res[1]) = get_zx_distances(mat, (COMPUTE_TYPE)compute_type);

    return res;
}

int get_distance(int n, int m, int compute_type, char **stab_mat) {
    return get_distance(
        from_c_array(n, m, stab_mat),
        (COMPUTE_TYPE)compute_type
    );
}
