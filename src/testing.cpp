#include <algorithm>
#include <fstream>
#include <iostream>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "brouwer_zimmerman.hpp"

BMatrix steane_code() {
    const int n = 6;
    const int m = 7;

    // 0 is I, 1 is Z, 2 is X, 3 is Y
    const int entries[n][m] = {
        {0, 0, 0, 1, 1, 1, 1},
        {0, 0, 0, 3, 3, 3, 3},
        {0, 1, 1, 0, 0, 1, 1},
        {0, 3, 3, 0, 0, 3, 3},
        {1, 0, 1, 0, 1, 0, 1},
        {3, 0, 3, 0, 3, 0, 3}
    };

    BMatrix steane(n, 2 * m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (entries[i][j] & 1)
                steane.set(i, 2 * j, 1);
            if (entries[i][j] & 2)
                steane.set(i, 2 * j + 1, 1);
        }
    }

    return steane;
}

BMatrix bmatrix_conversion(std::vector<std::string> code) {
    const int m = 2 * code.back().size();
    const int n = code.size();

    BMatrix result;

    for (int i = 0; i < n; ++i) {
        BVector row(m);

        for (int j = 0; j < code[i].size(); ++j) {
            switch (code[i][j]) {
                case '_': break;
                case 'I': break;
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

std::vector<BMatrix> ben_codes() {
    std::ifstream fi("../testing/d4_codes.txt");

    std::vector<BMatrix> bencodes;
    std::vector<std::string> code;
    std::string current;

    while (getline(fi, current)) {
        if (current == "") {
            bencodes.push_back(bmatrix_conversion(code));
            code.clear();
        }
        else {
            code.push_back(current);
        }
    }

    return bencodes;
}

int main() {
    auto codes = ben_codes();

    for (auto code: codes) {
        int z, x;
        std::tie(z, x) = get_zx_distances(code);
        std::cout << x << ' ' << z << std::endl;
    }

    return 0;
}
