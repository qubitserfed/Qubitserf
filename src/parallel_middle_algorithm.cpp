#include "parallel_middle_algorithm.hpp"

using u32 = unsigned int;

int parallel_middle_algorithm(BMatrix stab_mat, BMatrix anticomms, COMPUTE_TYPE compute_type, bool verbose_flag) {
    Printer printer(verbose_flag);
    const int m = stab_mat.m / 2;

    ParallelHashTable *s0, *s1;
    s0 = new ParallelHashTable();
    s1 = new ParallelHashTable();

    s0->insert(BVector(stab_mat.n), BVector(anticomms.n));

    for (int d = 1; d <= m; ++d) {
        std::swap(s0, s1);
        
        s0->reset_for_code(m, d, false);
        bool res = false;

        res = parallel_symplectic_combinations(m, d, [&] (BVector &v) -> bool {
            BVector stab_syndome, anticomms_syndrome;

            stab_syndome = BVector(stab_mat.n);
            for (int i = 0; i < stab_mat.n; ++i)
                stab_syndome.set(i, sym_prod(stab_mat.row(i), v));

            anticomms_syndrome = BVector(anticomms.n);
            for (int i = 0; i < anticomms.n; ++i)
                anticomms_syndrome.set(i, sym_prod(anticomms.row(i), v));

            // check for compatilble words of weight d-1
            if (s1->search(stab_syndome, anticomms_syndrome))
                return true;
            s0->insert(stab_syndome, anticomms_syndrome);
            return false;
        }, compute_type.no_threads);

        if (res) {
            delete s0;
            delete s1;
            printer("Distance: =", 2 * d - 1, "\nElapsed:", Timestamp(), "\n", ResetTime());
            return 2 * d - 1;
        }

        printer("Distance bound: >", 2 * d - 1, "\nElapsed:", Timestamp(), "\n", ResetTime());

        res = parallel_symplectic_combinations(m, d, [&] (BVector &v) -> bool {
            BVector stab_syndome, anticomms_syndrome;

            stab_syndome = BVector(stab_mat.n);
            for (int i = 0; i < stab_mat.n; ++i)
                stab_syndome.set(i, sym_prod(stab_mat.row(i), v));

            anticomms_syndrome = BVector(anticomms.n);
            for (int i = 0; i < anticomms.n; ++i)
                anticomms_syndrome.set(i, sym_prod(anticomms.row(i), v));

            if (s0->search(stab_syndome, anticomms_syndrome))
                return true;
            return false;
        }, compute_type.no_threads);

        if (res) {
            delete s0;
            delete s1;
            printer("Distance: =", 2 * d, "\nElapsed:", Timestamp(), "\n", ResetTime());
            return 2 * d;
        }

        printer("Distance bound: >", 2 * d, "\nElapsed:", Timestamp(), "\n", ResetTime());
    }

    delete s0;
    delete s1;
    return 1;
}

int parallel_css_middle_algorithm(BMatrix stab_mat, BMatrix anticomms, COMPUTE_TYPE compute_type, bool verbose_flag) {
    Printer printer(verbose_flag);
    const int m = stab_mat.m / 2;

    ParallelHashTable *s0, *s1;
    s0 = new ParallelHashTable();
    s1 = new ParallelHashTable();

    s0->insert(BVector(stab_mat.n), BVector(anticomms.n));

    for (int d = 1; d <= m; ++d) {
        std::swap(s0, s1);
        
        s0->reset_for_code(m, d, false);
        bool res = false;

        res = parallel_combinations(m, d, [&] (BVector &v) -> bool {
            BVector stab_syndome, anticomms_syndrome;

            stab_syndome = transposed_product(v, stab_mat);
            anticomms_syndrome = transposed_product(v, anticomms);

            // check for compatilble words of weight d-1
            if (s1->search(stab_syndome, anticomms_syndrome))
                return true;
            s0->insert(stab_syndome, anticomms_syndrome);
            return false;
        }, compute_type.no_threads);

        if (res) {
            delete s0;
            delete s1;
            printer("Distance: =", 2 * d - 1, "\nElapsed:", Timestamp(), "\n", ResetTime());
            return 2 * d - 1;
        }

        printer("Distance bound: >", 2 * d - 1, "\nElapsed:", Timestamp(), "\n", ResetTime());

        res = parallel_combinations(m, d, [&] (BVector &v) -> bool {
            BVector stab_syndome, anticomms_syndrome;

            stab_syndome = transposed_product(v, stab_mat);
            anticomms_syndrome = transposed_product(v, anticomms);

            if (s0->search(stab_syndome, anticomms_syndrome))
                return true;
            return false;
        }, compute_type.no_threads);

        if (res) {
            delete s0;
            delete s1;
            printer("Distance: =", 2 * d, "\nElapsed:", Timestamp(), "\n", ResetTime());
            return 2 * d;
        }

        printer("Distance bound: >", 2 * d, "\nElapsed:", Timestamp(), "\n", ResetTime());
    }

    delete s0;
    delete s1;
    return 1;
}
