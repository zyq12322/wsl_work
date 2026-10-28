#include "process_input_var.h"
#include "variables.h"
#include <iostream>

int main() {
    try {
        // 1. 读取输入文件
        process_input_var::get_IV("DNO_input.txt");

        // 2. 检查输入合法性
        process_input_var::check_IV();

        // 3. 如果启用无量纲化，执行无量纲化
        if (nondimensionalize_flag) {
            process_input_var::nondimensionalize();
        }

        // 4. 打印解析后的变量
        process_input_var::print_IV();

        // 5. 生成配置文件（用于验证输出）
        writeConfigFile("test_output_config.txt");
        std::cout << "\nConfig file written to test_output_config.txt\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
