#include "var_op_function.h"

#include <iostream>
//#include <fstream>

int main()
{
    using namespace var_opf;

    try
    {
        std::cout << "=============================\n";
        std::cout << " trim()\n";
        std::cout << "=============================\n";

        std::string s1 = "   Hello World   ";
        std::cout << "[" << trim(s1) << "]\n\n";


        std::cout << "=============================\n";
        std::cout << " to_lower()\n";
        std::cout << "=============================\n";

        std::cout << to_lower("HeLLo WoRLD") << "\n\n";


        std::cout << "=============================\n";
        std::cout << " is_integer()\n";
        std::cout << "=============================\n";

        std::cout << "123      : " << is_integer("123") << '\n';
        std::cout << "-456     : " << is_integer("-456") << '\n';
        std::cout << "+789     : " << is_integer("+789") << '\n';
        std::cout << "12.3     : " << is_integer("12.3") << '\n';
        std::cout << "abc      : " << is_integer("abc") << "\n\n";


        std::cout << "=============================\n";
        std::cout << " is_float()\n";
        std::cout << "=============================\n";

        std::cout << "123       : " << is_float("123") << '\n';
        std::cout << "3.14159   : " << is_float("3.14159") << '\n';
        std::cout << "-1.2e-3   : " << is_float("-1.2e-3") << '\n';
        std::cout << "1e10      : " << is_float("1e10") << '\n';
        std::cout << "abc       : " << is_float("abc") << "\n\n";


        std::cout << "=============================\n";
        std::cout << " parse_string_toIDBS()\n";
        std::cout << "=============================\n";

        IDBSVariant v1 = parse_string_toIDBS("123");
        IDBSVariant v2 = parse_string_toIDBS("3.1415926535");
        IDBSVariant v3 = parse_string_toIDBS("true");
        IDBSVariant v4 = parse_string_toIDBS("hello world");

        std::cout << variantIDBS_to_string(v1) << '\n';
        std::cout << variantIDBS_to_string(v2) << '\n';
        std::cout << variantIDBS_to_string(v3) << '\n';
        std::cout << variantIDBS_to_string(v4) << "\n\n";


        std::cout << "=============================\n";
        std::cout << " variantIDBS_to_double()\n";
        std::cout << "=============================\n";

        std::cout << variantIDBS_to_double(v1) << '\n';
        std::cout << variantIDBS_to_double(v2) << '\n';
        std::cout << variantIDBS_to_double(v3) << "\n\n";


        std::cout << "=============================\n";
        std::cout << " split_by_comma()\n";
        std::cout << "=============================\n";

        std::vector<std::string> tokens;

        split_by_comma(
            "apple, banana, orange , grape",
            tokens
        );

        for (const auto& token : tokens)
            std::cout << token << '\n';

        std::cout << '\n';


        std::cout << "=============================\n";
        std::cout << " get_key_before_char()\n";
        std::cout << " get_val_after_char()\n";
        std::cout << "=============================\n";

        std::string line = "temperature = 273.15";

        std::cout
            << "key   = "
            << get_key_before_char(line,'=')
            << '\n';

        std::cout
            << "value = "
            << get_val_after_char(line,'=')
            << "\n\n";


        std::cout << "=============================\n";
        std::cout << " print_variantIDBS_vector()\n";
        std::cout << "=============================\n";

        IDBSVector vec =
        {
            100,
            3.14,
            true,
            std::string("example")
        };

        print_variantIDBS_vector(
            std::cout,
            vec,
            " | "
        );

        std::cout << "\n\n";


        std::cout << "=============================\n";
        std::cout << " divide_at_indexofIDBS()\n";
        std::cout << "=============================\n";

        divide_at_indexofIDBS(vec, 2.0, 0);
        divide_at_indexofIDBS(vec, 2.0, 1);

        print_variantIDBS_vector(
            std::cout,
            vec,
            " | "
        );

        std::cout << "\n\n";


        std::cout << "=============================\n";
        std::cout << " Filesystem Test\n";
        std::cout << "=============================\n";

        check_or_create_output_directory(
            "test_output"
        );

        std::cout
            << "Directory check passed.\n\n";


        std::cout << "=============================\n";
        std::cout << " Exception Tests\n";
        std::cout << "=============================\n";

        try
        {
            divide_at_indexofIDBS(
                vec,
                0.0,
                0
            );
        }
        catch(const std::exception& e)
        {
            std::cout
                << "Expected exception: "
                << e.what()
                << '\n';
        }

        try
        {
            variantIDBS_to_double(
                std::string("abc")
            );
        }
        catch(const std::exception& e)
        {
            std::cout
                << "Expected exception: "
                << e.what()
                << '\n';
        }

        std::cout << "\nAll tests completed.\n";
    }
    catch(const std::exception& e)
    {
        std::cerr
            << "\nFatal error: "
            << e.what()
            << '\n';

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}