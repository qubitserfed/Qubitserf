#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <utility>

#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "middle_algorithm.hpp"

struct SyndromeCell {
    int population_count;
    std::map< BVector, BVector > anticomm_syndromes;

    SyndromeCell() {
        population_count = 0;
    }
};

int middle_algorithm(BMatrix stab_mat, BMatrix anticomms) {
    const int m = stab_mat.m;

    std::map< BVector, std::set<BVector> > s0, s1;
    s0[BVector(stab_mat.n)].insert(BVector(anticomms.n));

    for (int d = 1; d <= m; ++d) {
        std::swap(s0, s1);
        s0.clear();

        try {
            combinations(m, d, [&] (const std::vector<bool> &v0) {
                BVector stab_syndome, anticomms_syndrome, v;
                std::map< BVector, std::set<BVector> >::iterator it;

                v = v0;
                stab_syndome = transposed_product(v, stab_mat);
                anticomms_syndrome = transposed_product(v, anticomms);

                // check for compatilble words of weight d-1
                it = s1.find(stab_syndome);
                if (it != s1.end()) {
                    if (it->second.size() > 2) {
                        throw 2 * d - 1;
                    }
                    else if (*it->second.begin() != anticomms_syndrome) {
                        throw 2 * d - 1;
                    }
                }

                it = s0.find(stab_syndome);
                if (it != s0.end()) {
                    if (it->second.size() > 2) {
                        throw 2 * d;
                    }
                    else if (*it->second.begin() != anticomms_syndrome) {
                        throw 2 * d;
                    }
                    else {
                        it->second.insert(anticomms_syndrome);
                    }
                }
                else {
                    s0[stab_syndome].insert(anticomms_syndrome);
                }
            });
        }
        catch (int result) {
            return result;
        }
    }

    my_assert(0);
    return 0;
}
