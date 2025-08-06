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
            if (s1->find(stab_syndome, anticomms_syndrome))
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

            if (s0->find(stab_syndome, anticomms_syndrome))
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
    const int m = stab_mat.m;
    int l, r;

    l = m / 2;
    r = (m + 1) / 2;

    Printer printer(verbose_flag);

    ParallelHashTable *sl, *sr;
    sl = new ParallelHashTable();
    sr = new ParallelHashTable();

    sl->reset_for_code(l, 25, true); // implement a resize for code so that it can be resized after each iteration
    sr->reset_for_code(r, 25, true);

    sl->insert(BVector(stab_mat.n), BVector(anticomms.n));
    sr->insert(BVector(stab_mat.n), BVector(anticomms.n));

    for (int d = 1; d <= m; ++d) {
        bool found = false;

        // 1) search for a compatible word of weight 2*d-1
        // 1.i) search in left half
        found = parallel_combinations(l, d, [&] (BVector &v) {
            BVector w = v;
            w.resize(m); // w is in the left half, its pair is searched for in the right half
            BVector stab_syndome = transposed_product(w, stab_mat), anticomms_syndrome = transposed_product(w, anticomms);
            return sr->find(stab_syndome, anticomms_syndrome);
        }, compute_type.no_threads) ||
        parallel_combinations(r, d, [&] (BVector &v) { // 1.ii) search in right half
            BVector w = v;
            w.resize(m);
            w = shift_left(w, l); // w is in the right half, its pair is searched for in the left half
            BVector stab_syndome = transposed_product(w, stab_mat), anticomms_syndrome = transposed_product(w, anticomms);
            return sl->find(stab_syndome, anticomms_syndrome);
        }, compute_type.no_threads);

        if (found) {
            delete sl;
            delete sr;
            printer("Distance: =", 2 * d - 1, "\nElapsed:", Timestamp(), "\n", ResetTime());
            return 2 * d - 1;
        }

        printer("Distance bound: >", 2 * d - 1, "\nElapsed:", Timestamp(), "\n", ResetTime());

        // 2) populate halves with the weight d words
        // 2.i) add to left half
        parallel_combinations(l, d, [&] (BVector &v) {
            BVector w = v;
            w.resize(m);
            BVector stab_syndome = transposed_product(w, stab_mat), anticomms_syndrome = transposed_product(w, anticomms);
            sl->insert(stab_syndome, anticomms_syndrome);
            return false;
        }, compute_type.no_threads);

        // 2.ii) add to right half
        parallel_combinations(r, d, [&] (BVector &v) {
            BVector w = v;
            w.resize(m);
            w = shift_left(w, l);
            BVector stab_syndome = transposed_product(w, stab_mat), anticomms_syndrome = transposed_product(w, anticomms);
            sr->insert(stab_syndome, anticomms_syndrome);
            return false;
        }, compute_type.no_threads);


        // 3) search for a compatible word of weight 2*d
        found = parallel_combinations(l, d, [&] (BVector &v) { // 3.i) search in left half
            BVector w = v;
            w.resize(m);
            BVector stab_syndome = transposed_product(w, stab_mat), anticomms_syndrome = transposed_product(w, anticomms);
            return sr->find(stab_syndome, anticomms_syndrome);
        }, compute_type.no_threads) ||
        parallel_combinations(r, d, [&] (BVector &v) { // 3.ii) search in right half
            BVector w = v;
            w.resize(m);
            w = shift_left(w, l); // w is in the right half, its pair is searched for in the left half
            BVector stab_syndome = transposed_product(w, stab_mat), anticomms_syndrome = transposed_product(w, anticomms);
            return sl->find(stab_syndome, anticomms_syndrome);
        }, compute_type.no_threads);

        if (found) {
            delete sl;
            delete sr;
            printer("Distance: =", 2 * d, "\nElapsed:", Timestamp(), "\n", ResetTime());
            return 2 * d;
        }

        printer("Distance bound: >", 2 * d, "\nElapsed:", Timestamp(), "\n", ResetTime());
    }

    delete sl;
    delete sr;
    return 1;
}