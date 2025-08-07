#include "operator_weight.hpp"

#include "linear_algebra.hpp"
#include "quantum_utilities.hpp"
#include "combinatorics.hpp"
#include "utility.hpp"

int get_css_operator_weight(BMatrix stabs, BMatrix logs, BMatrix destabs, BVector vec, COMPUTE_TYPE compute_type = {true, false, 1}, bool verbose = false) {
    const int n = stabs.m;
    Printer printer(verbose);
    ParallelHashTable *s0, *s1;

    if (stabs.n == 0) {
        printer("Distance: =", vec.weight(), "\n");
        return vec.weight();
    }
    if (in_span(stabs, vec)) {
        printer("Distance: =0\n");
        return 0;
    }

    s0 = new ParallelHashTable();
    s1 = new ParallelHashTable();

    BVector vec_stab_syndrome = transposed_product(vec, stabs);
    BVector vec_log_syndrome = transposed_product(vec, logs);

    printer(ResetTime());

    s0->insert(BVector(logs.n), BVector(logs.m));
    for (int d = 1; d <= n; ++d) {
        std::swap(s0, s1);
        s0->reset_for_code(n, d, false);

        bool found = false;    

        found = parallel_combinations(n, d, [&](BVector cand) -> bool{
            BVector cand_stab_syndrome = transposed_product(cand, stabs);
            BVector cand_log_syndrome = transposed_product(cand, logs);

            s0->insert(cand_log_syndrome, cand_stab_syndrome);
            return s1->find(cand_log_syndrome + vec_log_syndrome, cand_stab_syndrome + vec_stab_syndrome);
        }, compute_type.no_threads);

        if (found) {
            printer("Distance: =", 2 * d - 1, "\nElapsed:", Timestamp(), "\n", ResetTime());
            delete s0;
            delete s1;
            return 2 * d - 1;
        }
        printer("Distance bound: >", 2 * d - 1, "\nElapsed:", Timestamp(), "\n", ResetTime());

        found = parallel_combinations(n, d, [&](BVector cand) -> bool{
            BVector cand_stab_syndrome = transposed_product(cand, stabs);
            BVector cand_log_syndrome = transposed_product(cand, logs);

            return s0->find(cand_log_syndrome + vec_log_syndrome, cand_stab_syndrome + vec_stab_syndrome);
        }, compute_type.no_threads);

        if (found) {
            printer("Distance: =", 2 * d, "\nElapsed:", Timestamp(), "\n", ResetTime());
            delete s0;
            delete s1;
            return 2 * d;
        }
        printer("Distance bound: >", 2 * d, "\nElapsed:", Timestamp(), "\n", ResetTime());
    }

//    my_assert(false); // should never happen

    delete s0;
    delete s1;
    return 0;
}

std::pair<int, int> get_zx_operator_weight(BMatrix stabs, BVector vec, COMPUTE_TYPE compute_type = {true, false, 1}, bool verbose = false) {
    Printer printer(verbose);
    BMatrix logs = logical_operators(stabs);
    BMatrix destabs = destabilizers(stabs, logs);

    BVector z_vec, x_vec;
    std::tie(z_vec, x_vec) = zx_parts(vec);

    BMatrix z_ops, x_ops;
    std::tie(z_ops, x_ops) = zx_parts(logs);
    z_ops.remove_zeros(); x_ops.remove_zeros();

    BMatrix z_stabs, x_stabs;
    std::tie(z_stabs, x_stabs) = zx_parts(stabs);
    z_stabs.remove_zeros(); x_stabs.remove_zeros();

    BMatrix z_destabs, x_destabs;
    std::tie(z_destabs, x_destabs) = zx_parts(destabs);
    z_destabs.remove_zeros(); x_destabs.remove_zeros();

    printer("Computing Z-weight modulo stabilizers...\n");
    int z_dist = get_css_operator_weight(z_stabs, x_ops, x_destabs, z_vec, compute_type, verbose);
    printer("Computing X-weight modulo stabilizers...\n");
    int x_dist = get_css_operator_weight(x_stabs, z_ops, z_destabs, x_vec, compute_type, verbose);
    printer("Total Elapsed:", Timestamp(), "\n");

    return std::make_pair(z_dist, x_dist);
}

int get_operator_weight(BMatrix stabs, BVector vec, COMPUTE_TYPE compute_type = {true, false, 1}, bool verbose = false) {
    if (is_css(stabs)) {
        int x_dist, z_dist;
        std::tie(x_dist, z_dist) = get_zx_operator_weight(stabs, vec, compute_type, verbose);
        return std::max(x_dist, z_dist);
    }
    else {
        std::cerr << "Operator weight not available for non-css codes" << std::endl;
        my_assert(0);
    }
}
