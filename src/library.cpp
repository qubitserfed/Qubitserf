#include <tuple>
#include <utility>

#include "linear_algebra.hpp"
#include "brouwer_zimmerman.hpp"

extern "C" {
    int *get_zx_distances(int, int, int, char **);
    int get_distance(int, int, int, char **);
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
