#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <cstring>

#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "quantum_utilities.hpp"
#include "brouwer_zimmerman.hpp"

using u64 = unsigned long long;

BMatrix bmatrix_conversion(std::vector<std::string> code) {
    const int m = 2 * code.back().size();
    const int n = code.size();

    BMatrix result;

    for (int i = 0; i < n; ++i) {
        BVector row(m);

        for (int j = 0; j < code[i].size(); ++j) {
            switch (code[i][j]) {
                case '.':
                case '_':
                case 'I':
                    break;
                case 'X':
                    row.set(2 * j + 1, 1);
                    break;
                case 'Z':
                    row.set(2 * j, 1);
                    break;
                case 'Y':
                    row.set(2 * j, 1);
                    row.set(2 * j + 1, 1);
                    break;
                default:
                    my_assert(0);
                    break;
            }
        }

        result.append_row(row);
    }

    return result;
}

int main(int argc, char **argv) {
    std::string current;
    std::vector<std::string> code_strs;
    
    while (getline(std::cin, current)) {
        if (current == "")
            break;
        code_strs.push_back(current);
    }

    BMatrix code = bmatrix_conversion(code_strs);

    if (argc > 0 && !strcmp(argv[1], "zx")) {
        int z_dist, x_dist;
        
        std::tie(z_dist, x_dist) = get_zx_distances(code);
        std::cout << z_dist << ' ' << x_dist << std::endl;
    }
    else {
        std::cout << get_distance(code) << std::endl;
    }
    return 0;
}
