#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <cstring>

#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "quantum_utilities.hpp"
#include "brouwer_zimmerman.hpp"
#include "middle_algorithm.hpp"

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
    std::vector< std::string > argv_flags;
    COMPUTE_TYPE compute_type;
    std::vector<std::string>::iterator it;
    bool zx_flag = false, middle_flag = false;
    
    compute_type = (COMPUTE_TYPE) {
        true,
        false,
        1
    };

    while (getline(std::cin, current)) {
        if (current == "")
            break;
        code_strs.push_back(current);
    }

    BMatrix code = bmatrix_conversion(code_strs);

    for (int i = 1; i < argc; ++i)
        argv_flags.push_back( std::string(argv[i]) );


    it = std::find(argv_flags.begin(), argv_flags.end(), "--middle");
    if (it != argv_flags.end())
        middle_flag = true;

    it = std::find(argv_flags.begin(), argv_flags.end(), "--threads");
    if (it != argv_flags.end()) {
        it++;
        if (it != argv_flags.end()) {
            compute_type.no_threads = std::stoi(*it);
        }
        else {
            std::cerr << "--threads flag need be followed by the number of threads!\n";
        }
    }

    if (std::find(argv_flags.begin(), argv_flags.end(), "--zx") != argv_flags.end())
        zx_flag = true;

    if (zx_flag) {
        int z_dist, x_dist;
        
        if (middle_flag)
            std::tie(z_dist, x_dist) = get_zx_distances_with_middle(code);
        else
            std::tie(z_dist, x_dist) = get_zx_distances(code, compute_type);
        std::cout << z_dist << ' ' << x_dist << std::endl;
    }
    else {
        if (middle_flag)
            std::cout << get_distance_with_middle(code) << std :: endl;
        else
            std::cout << get_distance(code, compute_type) << std::endl;
    }

    return 0;
}
