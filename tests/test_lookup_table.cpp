#include "catch.hpp"
#include "../src/table_library.hpp"
#include <iostream>

TEST_CASE("Steane test") {
    char* inputs[] = {
        (char*)"XXXXIII",
        (char*)"IXXIXXI",
        (char*)"IIXXIXX",
        (char*)"ZZZZIII",
        (char*)"IZZIZZI",
        (char*)"IIZZIZZ"
    };

    int n = 7;
    int k = 1;

    CSS_Code* code = css_parse_code(inputs, n, k);

    std::cout << "Z Stabilizers:" << std::endl;
    print(code->z_stabs);
    std::cout << "X Stabilizers:" << std::endl;
    print(code->x_stabs);

    CSS_Table *table = css_make_table(n, k, 3, code);

    auto word = parse_word("XXIIXII", n);

    auto res = css_lookup_word(word, code, table);
    
    std::cout << res.first << " " << res.second <<std::endl;

    // Cleanup
    delete code;
    delete table;
}

/*
stabs in row-echelon form:
XIIXXXI
IXIXXIX
IIXXIXX

random error:
IIIXIXI

table->x_lookup[2][1010] == true if and only if the error with reduced form 1010 is equivalent to a weight 2 error 

X stabcomms + X logicals + X anticomms
Z anticomms + Z logicals + Z stabs
*/
