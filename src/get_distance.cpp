#include "get_distance.hpp"

#include "singlethread_middle_algorithm.hpp"
#include "parallel_middle_algorithm.hpp"
#include "brouwer_zimmerman.hpp"

#include <iostream>
#include <cstdlib>


int get_distance(BMatrix stab_mat, ALGORITHMS algorithm, COMPUTE_TYPE compute_type, bool verbose_flag) {
    if (is_css(stab_mat)) {
        int x_dist, z_dist;
        std::tie(x_dist, z_dist) = get_zx_distances(stab_mat, algorithm, compute_type, verbose_flag);
        return std::min(x_dist, z_dist);
    }
    else {
        if (algorithm == MIDDLE_ALGORITHM) {
            if (compute_type.no_threads == 1)
                return middle_algorithm(stab_mat, logical_operators(stab_mat), verbose_flag);
            else
                return parallel_middle_algorithm(stab_mat, logical_operators(stab_mat), compute_type, verbose_flag);
        }
        else if (algorithm == BROUWER_ZIMMERMAN_ALGORITHM) {
            std::cerr << "Brouwer-Zimmerman algorithm is not implemented for non-CSS codes.\n";
            exit(1);
        }
    }
}

std::pair<int, int> get_zx_distances(BMatrix stab_mat, ALGORITHMS algorithm, COMPUTE_TYPE compute_type, bool verbose_flag) {
    BMatrix closure_mat, x_stab, z_stab, x_ops, z_ops;
    Printer printer(verbose_flag);

    closure_mat = logical_operators(stab_mat);

    std::tie(z_ops, x_ops) = zx_parts(closure_mat);
    std::tie(z_stab, x_stab) = zx_parts(stab_mat);

    to_row_echelon(z_stab);
    to_row_echelon(x_stab);

    z_stab.remove_zeros();
    x_stab.remove_zeros();

    z_ops.remove_zeros();
    x_ops.remove_zeros();

    if (algorithm == MIDDLE_ALGORITHM) {
        if (compute_type.no_threads == 1) {
            printer("Z-distance:\n");
            const int z_dist = css_middle_algorithm(z_stab, z_ops, verbose_flag);
            printer("X-distance:\n");
            const int x_dist = css_middle_algorithm(x_stab, x_ops, verbose_flag);            
            return std::make_pair(z_dist, x_dist);
        }
        else {
            printer("Z-distance:\n");
            const int z_dist = parallel_css_middle_algorithm(z_stab, z_ops, compute_type, verbose_flag);
            printer("X-distance:\n");
            const int x_dist = parallel_css_middle_algorithm(x_stab, x_ops, compute_type, verbose_flag);
            return std::make_pair(z_dist, x_dist);
        }
    }
    else {
        printer("Z-distance:\n");
        const int z_dist = brouwer_zimmerman(z_stab, z_ops, compute_type, verbose_flag);
        printer("X-distance:\n");
        const int x_dist = brouwer_zimmerman(x_stab, x_ops, compute_type, verbose_flag);
        return std::make_pair(z_dist, x_dist);
    }
}
