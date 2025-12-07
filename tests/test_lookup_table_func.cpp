#include "catch.hpp"
#include "../src/table_library.hpp"
#include <iostream>

TEST_CASE("Steane test") {
    // [IIIXXXX IXXIXXI XIXIXIX IIIZZZZ IZZIZZI ZIZIZIZ]
    char* inputs[] = {
        (char*)"XXXXIII",
        (char*)"IXXIXXI",
        (char*)"XIXIXIX",
        (char*)"IIIZZZZ",
        (char*)"IZZIZZI",
        (char*)"ZIZIZIZ"
    };

    int n = 7;
    int k = 1;

    CSS_Code* code = css_parse_code(inputs, n, k);

    std::cout << "Z Stabilizers:" << std::endl;
    print(code->z_stabs);
    std::cout << "X Stabilizers:" << std::endl;
    print(code->x_stabs);

    CSS_Table *table = css_make_table(n, k, 3, code);

    std::cout << "parsing word" << std::endl;
    auto word = parse_word("XXXIIII", n);
    print(word);

    std::cout << "lookup word" << std::endl;
    auto res = css_lookup_word(word, code, table);
    std::cout<<res.first<<" "<<res.second<<std::endl;

    // Cleanup
    delete code;
    delete table;
}
