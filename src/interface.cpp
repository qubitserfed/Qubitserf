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
    bool zx_flag = false, bz_flag = false, verbose_flag = false, parallel_flag = false;

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

    it = std::find(argv_flags.begin(), argv_flags.end(), "--bz");
    if (it != argv_flags.end()) {
        if (!is_css(code)) {
            std :: cerr << "The Brouwer-Zimmerman Algorithm does not currently support non-CSS codes!\n";
            exit(1);
        }
        bz_flag = true;
    }

    it = std::find(argv_flags.begin(), argv_flags.end(), "--zx");
    if (it != argv_flags.end()) {
        if (!is_css(code)) {
            std :: cerr << "You cannot ask for the Z and X distances of a non-CSS code.\n";
            exit(1);
        }
        zx_flag = true;
    }

    it = std::find(argv_flags.begin(), argv_flags.end(), "--threads");
    if (it != argv_flags.end()) {
        parallel_flag = true;

        it++;
        if (it != argv_flags.end()) {
            // Check that *it is a number
            if ((*it).empty() || !std::all_of((*it).begin(), (*it).end(), ::isdigit) || std::stoi(*it) <= 0) {
                std::cerr << "--threads flag must be followed by a positive integer!\n";
                exit(1);
            }
            compute_type.no_threads = std::stoi(*it);
        }
        else {
            std::cerr << "--threads flag need be followed by the number of threads!\n";
            exit(1);
        }
    }

    it = std::find(argv_flags.begin(), argv_flags.end(), "--verbose");
    if (it != argv_flags.end()) {
        verbose_flag = true;
    }

    it = std::find(argv_flags.begin(), argv_flags.end(), "-v");
    if (it != argv_flags.end()) {
        verbose_flag = true;
    }

    if (zx_flag) {
        int z_dist, x_dist;
        if (parallel_flag) {
            if (bz_flag) {
                if (verbose_flag)
                    get_zx_distances(code, compute_type, verbose_flag);
                else {
                    std::tie(z_dist, x_dist) = get_zx_distances(code, compute_type, verbose_flag);
                    std::cout << z_dist << ' ' << x_dist << std::endl;
                }
            }
            else {
                if (verbose_flag)
                    get_zx_distances_with_parallelized_middle(code, compute_type, verbose_flag);
                else {
                    std::tie(z_dist, x_dist) = get_zx_distances_with_parallelized_middle(code, compute_type, verbose_flag);
                    std::cout << z_dist << ' ' << x_dist << std::endl;
                }
            }
        }
        else {
            if (bz_flag) {
                if (verbose_flag)
                    get_zx_distances(code, compute_type, verbose_flag);
                else {
                    std::tie(z_dist, x_dist) = get_zx_distances(code, compute_type, verbose_flag);
                    std::cout << z_dist << ' ' << x_dist << std::endl;
                }
            }
            else {
                if (verbose_flag)
                    get_zx_distances_with_middle(code, verbose_flag);
                else {
                    std::tie(z_dist, x_dist) = get_zx_distances_with_middle(code, verbose_flag);
                    std::cout << z_dist << ' ' << x_dist << std::endl;
                }
            }
        }
    }
    else {
        if (parallel_flag) {
            if (bz_flag) {
                if (verbose_flag)
                    get_distance(code, compute_type, verbose_flag);
                else
                    std::cout << get_distance(code, compute_type, verbose_flag) << std::endl;
            }
            else {
                if (verbose_flag)
                    get_distance_with_parallelized_middle(code, compute_type, verbose_flag);
                else
                    std::cout << get_distance_with_parallelized_middle(code, compute_type, verbose_flag) << std :: endl;
            }
        }
        else {
            if (bz_flag) {
                if (verbose_flag)
                    get_distance(code, compute_type, verbose_flag);
                else
                    std::cout << get_distance(code, compute_type, verbose_flag) << std::endl;
            }
            else {
                if (verbose_flag)
                    get_distance_with_middle(code, verbose_flag);
                else
                    std::cout << get_distance_with_middle(code, verbose_flag) << std :: endl;
            }
        }
    }

    return 0;
}
