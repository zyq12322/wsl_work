#include "process_input_var.h"
#include "var_op_function.h"
#include "variables.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <muParser.h>



namespace process_input_var {

    void get_GIV(const std::string& trimmed_line){
        // input the trimmed line string with format "var_key = var_value".
        std::string key = var_opf::get_key_before_char(trimmed_line, '=');

        if (key == var_opf::trim(trimmed_line)) {
            std::cerr << "Warning: Skipping line without '=': " << trimmed_line << std::endl;
            return;
        }

        if (key.empty()) {
            std::cerr << "Warning: Skipping line without var_key: " << trimmed_line << std::endl;
            return;
        }
        
        // return trimmed var_value string.
        std::string value = var_opf::get_val_after_char(trimmed_line, '=');

        if (key == "dimension_num") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping GIV.dimension_num = " << GIV.dimension_num << std::endl;
                return;
            }
            int temp;
            if (var_opf::try_parse_int(value, temp)) {
                GIV.dimension_num = temp;
            } else {
                std::cerr << "Error: Integer out of range or conversion failed for '" << value 
                << "', using default GIV.dimension_num = " << GIV.dimension_num << std::endl;
            }
            return;

        }
        else if (key == "DNO_M") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping GIV.DNO_M = " << GIV.DNO_M << std::endl;
                return;
            }
            int temp;
            if (var_opf::try_parse_int(value, temp)) {
                GIV.DNO_M = temp;
            } else {
                std::cerr << "Error: Integer out of range or conversion failed for '" << value 
                << "', using default GIV.DNO_M = " << GIV.DNO_M << std::endl;
            }
            return;
            
        }
        else if (key == "DNO_Mb") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping GIV.DNO_Mb = " << GIV.DNO_Mb << std::endl;
                return;
            }
            int temp;
            if (var_opf::try_parse_int(value, temp)) {
                GIV.DNO_Mb = temp;
            } else {
                std::cerr << "Error: Integer out of range or conversion failed for '" << value 
                << "', using default GIV.DNO_Mb = " << GIV.DNO_Mb << std::endl;
            }
            return;

        }
        else if (key == "Dealiase_flag") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping GIV.Dealiase_flag = " << GIV.Dealiase_flag << std::endl;
                return;
            }
            std::string lower_val = var_opf::to_lower(value);
            if (lower_val == "true") {
                GIV.Dealiase_flag = true;
            }
            else if (lower_val == "false") {
                GIV.Dealiase_flag = false;
            }
            else{
                std::cerr << 
                "Warning: Parsing var_value failed! Use default GIV.Dealiase_flag = "
                << GIV.Dealiase_flag << std::endl;
            }
            return;
        }
        else if (key == "t_total") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping GIV.t_total = " << GIV.t_total << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                GIV.t_total = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default GIV.t_total = " << GIV.t_total << std::endl;
            }
            return;
            
        }
        else if (key == "dt") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping GIV.dt = " << GIV.dt << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                GIV.dt = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default GIV.dt = " << GIV.dt << std::endl;
            }
            return;

        }
        else if (key == "mean_depth") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping GIV.mean_depth = " << GIV.mean_depth << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                GIV.mean_depth = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default GIV.mean_depth = " << GIV.mean_depth << std::endl;
            }
            return;

        }
        else if (key == "gravity") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping GIV.gravity = " << GIV.gravity << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                GIV.gravity = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default GIV.gravity = " << GIV.gravity << std::endl;
            }
            return;

        }
        else if (key == "filter") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << std::endl;
                std::cerr << "  keeping GIV.filter = ";
                var_opf::print_variantIDBS_vector(std::cerr, GIV.filter);
                std::cerr << std::endl;
                return;
            }

            std::vector<std::string> parts;
            var_opf::split_by_comma(value, parts);

            GIV.filter.clear();
            for (const auto& part : parts) {
                std::variant<int, double, bool, std::string> temp = var_opf::parse_string_toIDBS(part);
                if (std::holds_alternative<std::string>(temp)) {
                    GIV.filter.push_back(var_opf::to_lower(std::get<std::string>(temp)));
                }
                else {
                    GIV.filter.push_back(temp);
                }
            }
            return;
        }
        else if (key == "NL_Ta") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping GIV.NL_Ta = " << GIV.NL_Ta << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                GIV.NL_Ta = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default GIV.NL_Ta = " << GIV.NL_Ta << std::endl;
            }
            return;

        }
        else {
            std::cerr << "Error: Illegal var_key " << key << " in var_class global_input_var" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    void get_IV2D(const std::string& trimmed_line) {
        // input the trimmed line string with format "var_key = var_value".
        std::string key = var_opf::get_key_before_char(trimmed_line, '=');

        if (key == var_opf::trim(trimmed_line)) {
            std::cerr << "Warning: Skipping line without '=': " << trimmed_line << std::endl;
            return;
        }

        if (key.empty()) {
            std::cerr << "Warning: Skipping line without var_key: " << trimmed_line << std::endl;
            return;
        }
        

        std::string value = var_opf::get_val_after_char(trimmed_line, '=');

        if (key == "Lx") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping IV2D.Lx = " << IV2D.Lx << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                IV2D.Lx = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default IV2D.Lx = " << IV2D.Lx << std::endl;
            }
            return;

        }
        else if (key == "Nx") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping IV2D.Nx = " << IV2D.Nx << std::endl;
                return;
            }
            int temp;
            if (var_opf::try_parse_int(value, temp)) {
                IV2D.Nx = temp;
            } else {
                std::cerr << "Error: Integer out of range or conversion failed for '" << value 
                << "', using default IV2D.Nx = " << IV2D.Nx << std::endl;
            }
            return;

        }
        else if (key == "x0") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping IV2D.x0 = " << IV2D.x0 << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                IV2D.x0 = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default IV2D.x0 = " << IV2D.x0 << std::endl;
            }
            return;

        }
        else if (key == "initial_conditions") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << std::endl;
                std::cerr << "  keeping IV2D.initial_conditions = ";
                var_opf::print_variantIDBS_vector(std::cerr, IV2D.initial_conditions);
                std::cerr << std::endl;
                return;
            }

            std::vector<std::string> parts;
            var_opf::split_by_comma(value, parts);

            IV2D.initial_conditions.clear();
            for (const auto& part : parts) {
                std::variant<int, double, bool, std::string> temp = var_opf::parse_string_toIDBS(part);
                if (std::holds_alternative<std::string>(temp)) {
                    IV2D.initial_conditions.push_back(var_opf::to_lower(std::get<std::string>(temp)));
                }
                else {
                    IV2D.initial_conditions.push_back(temp);
                }
            }
            return;
        }
        else if (key == "boundary_conditions") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << std::endl;
                std::cerr << "  keeping IV2D.boundary_conditions = " << IV2D.boundary_conditions << std::endl;
                return;
            }

            IV2D.boundary_conditions = var_opf::to_lower(value);
            return;
        }
        else {
            std::cerr << "Error: Illegal var_key " << key << " in var_class input_var_2d" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    void get_IV3D(const std::string& trimmed_line) {
        // input the trimmed line string with format "var_key = var_value".
        std::string key = var_opf::get_key_before_char(trimmed_line, '=');

        if (key == var_opf::trim(trimmed_line)) {
            std::cerr << "Warning: Skipping line without '=': " << trimmed_line << std::endl;
            return;
        }

        if (key.empty()) {
            std::cerr << "Warning: Skipping line without var_key: " << trimmed_line << std::endl;
            return;
        }
        

        std::string value = var_opf::get_val_after_char(trimmed_line, '=');

        if (key == "Lx") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping IV3D.Lx = " << IV3D.Lx << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                IV3D.Lx = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default IV3D.Lx = " << IV3D.Lx << std::endl;
            }
            return;

        }
        else if (key == "Ly") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping IV3D.Ly = " << IV3D.Ly << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                IV3D.Ly = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default IV3D.Ly = " << IV3D.Ly << std::endl;
            }
            return;

        }
        else if (key == "Nx") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping IV3D.Nx = " << IV3D.Nx << std::endl;
                return;
            }
            int temp;
            if (var_opf::try_parse_int(value, temp)) {
                IV3D.Nx = temp;
            } else {
                std::cerr << "Error: Integer out of range or conversion failed for '" << value 
                << "', using default IV3D.Nx = " << IV3D.Nx << std::endl;
            }
            return;

        }
        else if (key == "Ny") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping IV3D.Ny = " << IV3D.Ny << std::endl;
                return;
            }
            int temp;
            if (var_opf::try_parse_int(value, temp)) {
                IV3D.Ny = temp;
            } else {
                std::cerr << "Error: Integer out of range or conversion failed for '" << value 
                << "', using default IV3D.Ny = " << IV3D.Ny << std::endl;
            }
            return;

        }
        else if (key == "x0") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping IV3D.x0 = " << IV3D.x0 << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                IV3D.x0 = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default IV3D.x0 = " << IV3D.x0 << std::endl;
            }
            return;

        }
        else if (key == "y0") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping IV3D.y0 = " << IV3D.y0 << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                IV3D.y0 = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default IV3D.y0 = " << IV3D.y0 << std::endl;
            }
            return;

        }
        else if (key == "initial_conditions") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << std::endl;
                std::cerr << "  keeping IV3D.initial_conditions = ";
                var_opf::print_variantIDBS_vector(std::cerr, IV3D.initial_conditions);
                std::cerr << std::endl;
                return;
            }

            std::vector<std::string> parts;
            var_opf::split_by_comma(value, parts);

            IV3D.initial_conditions.clear();
            for (const auto& part : parts) {
                std::variant<int, double, bool, std::string> temp = var_opf::parse_string_toIDBS(part);
                if (std::holds_alternative<std::string>(temp)) {
                    IV3D.initial_conditions.push_back(var_opf::to_lower(std::get<std::string>(temp)));
                }
                else {
                    IV3D.initial_conditions.push_back(temp);
                }
            }
            return;
        }
        else if (key == "boundary_conditions") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << std::endl;
                std::cerr << "  keeping IV3D.boundary_conditions = " << IV3D.boundary_conditions << std::endl;
                return;
            }

            IV3D.boundary_conditions = var_opf::to_lower(value);
            return;
        }
        else {
            std::cerr << "Error: Illegal var_key " << key << " in var_class input_var_3d" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    void get_IDCV(const std::string& trimmed_line) {
        // input the trimmed line string with format "var_key = var_value".
        std::string key = var_opf::get_key_before_char(trimmed_line, '=');

        if (key == var_opf::trim(trimmed_line)) {
            std::cerr << "Warning: Skipping line without '=': " << trimmed_line << std::endl;
            return;
        }

        if (key.empty()) {
            std::cerr << "Warning: Skipping line without var_key: " << trimmed_line << std::endl;
            return;
        }
        

        std::string value = var_opf::get_val_after_char(trimmed_line, '=');

        if (key == "L_char") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping IDCV.L_char = " << IDCV.L_char << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                IDCV.L_char = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default IDCV.L_char = " << IDCV.L_char << std::endl;
            }
            return;

        }
        else if (key == "T_char") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping IDCV.T_char = " << IDCV.T_char << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                IDCV.T_char = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default IDCV.T_char = " << IDCV.T_char << std::endl;
            }
            return;

        }
        else if (key == "M_char") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping IDCV.M_char = " << IDCV.M_char << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                IDCV.M_char = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default IDCV.M_char = " << IDCV.M_char << std::endl;
            }
            return;

        }
        else {
            std::cerr << "Error: Illegal var_key " << key << " in var_class input_dimensionless_char_var" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    void get_VB(const std::string& trimmed_line) {
        // input the trimmed line string with format "var_key = var_value".
        std::string key = var_opf::get_key_before_char(trimmed_line, '=');

        if (key == var_opf::trim(trimmed_line)) {
            std::cerr << "Warning: Skipping line without '=': " << trimmed_line << std::endl;
            return;
        }

        if (key.empty()) {
            std::cerr << "Warning: Skipping line without var_key: " << trimmed_line << std::endl;
            return;
        }
        

        std::string value = var_opf::get_val_after_char(trimmed_line, '=');

        if (key == "time_dependency") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping VB.time_dependency = " << VB.time_dependency << std::endl;
                return;
            }
            std::string lower_val = var_opf::to_lower(value);
            if (lower_val == "true") {
                VB.time_dependency = true;
            }
            else if (lower_val == "false") {
                VB.time_dependency = false;
            }
            else{
                std::cerr << "Warning: Parsing var_value failed! Use default VB.time_dependency = " << VB.time_dependency << std::endl;
            }
            return;

        }
        else if (key == "beta") {

            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << std::endl;
                std::cerr << "  keeping VB.beta = " << VB.beta << std::endl;
                return;
            }

            VB.beta = value;

            return;
        }
        else if (key == "gamma") {

            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << std::endl;
                std::cerr << "  keeping VB.gamma = " << VB.gamma << std::endl;
                return;
            }

            VB.gamma = value;

            return;
        }
        else {
            std::cerr << "Error: Illegal var_key " << key << " in var_class varying_bottom" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    void get_SDC(const std::string& trimmed_line) {
        // input the trimmed line string with format "var_key = var_value".
        std::string key = var_opf::get_key_before_char(trimmed_line, '=');

        if (key == var_opf::trim(trimmed_line)) {
            std::cerr << "Warning: Skipping line without '=': " << trimmed_line << std::endl;
            return;
        }

        if (key.empty()) {
            std::cerr << "Warning: Skipping line without var_key: " << trimmed_line << std::endl;
            return;
        }
        

        std::string value = var_opf::get_val_after_char(trimmed_line, '=');

        if (key == "f_out") {
            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << 
                "\n  keeping SDC.f_out = " << SDC.f_out << std::endl;
                return;
            }
            double temp;
            if (var_opf::try_parse_double(value, temp)) {
                SDC.f_out = temp;
            } else {
                std::cerr << "Error: Double out of range or conversion failed for '" << value 
                << "', using default SDC.f_out = " << SDC.f_out << std::endl;
            }
            return;

        }
        else if (key == "out_file_path") {

            if (value.empty()) {
                std::cerr << "Warning: Skipping line without var_value: " << trimmed_line << std::endl;
                std::cerr << "  keeping SDC.out_file_path = " << SDC.out_file_path << std::endl;
                return;
            }

            SDC.out_file_path = value;

            return;
        }
        else {
            std::cerr << "Error: Illegal var_key " << key << " in var_class save_data_control." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    void get_IV(const std::string& input_file) {
        // output information when this function are running.
        std::cout << "\nSTART get_IV()!\n" << std::endl;
        // input file stream
        std::ifstream file(input_file);
        // if cann't open the specified file, output an error message and exit.
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open file " << input_file << std::endl;
            std::exit(EXIT_FAILURE);
        }
        // record line number
        int line_number = 0;
        std::string line;
        while (std::getline(file, line)) {
            // line_number starts from one.
            ++line_number;
            std::string trimmed_line = var_opf::trim(var_opf::get_key_before_char(line, '#'));
            if (trimmed_line.empty()) {
                continue;
            }
            std::string key = var_opf::get_key_before_char(trimmed_line, ':');
            
            if (key.empty() || key == trimmed_line) {
                std::cerr << "Error: Illegal var_class name at line " << line_number << ": " << trimmed_line << std::endl;
                std::exit(EXIT_FAILURE);
            }

            if (!var_opf::get_val_after_char(trimmed_line, ':').empty()){
                std::cerr << "Error: var_class must not have anything after colon " << line_number << ": " << trimmed_line << std::endl;
                std::exit(EXIT_FAILURE);
            }

            if (key == "global_input_var" ){
                std::cout << "\nProcessing var_class:" << key << std::endl;
                while (true) {
                    std::streampos pos = file.tellg();
                    if (!std::getline(file, line)) {
                        break;
                    }
                    ++line_number;
                    trimmed_line = var_opf::trim(var_opf::get_key_before_char(line, '#'));

                    if (trimmed_line.empty()) {
                        continue;
                    }
                    key = var_opf::get_key_before_char(trimmed_line, ':');
                    if (key != trimmed_line) {
                        if (key.empty()) {
                            std::cerr << "Error: Illegal var_class name at line " << line_number << ": " << trimmed_line << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        if (!file.seekg(pos)) {
                            std::cerr << "Error: Failed to seek back in file" << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        --line_number;
                        break;
                    }

                    get_GIV(trimmed_line);
                }
            }
            else if (key == "input_var_2d") {
                std::cout << "\nProcessing var_class:" << key << std::endl;
                while (true) {
                    std::streampos pos = file.tellg();
                    if (!std::getline(file, line)) {
                        break;
                    }
                    ++line_number;
                    trimmed_line = var_opf::trim(var_opf::get_key_before_char(line, '#'));

                    if (trimmed_line.empty()) {
                        continue;
                    }
                    key = var_opf::get_key_before_char(trimmed_line, ':');
                    if (key != trimmed_line) {
                        if (key.empty()) {
                            std::cerr << "Error: Illegal var_class name at line " << line_number << ": " << trimmed_line << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        if (!file.seekg(pos)) {
                            std::cerr << "Error: Failed to seek back in file" << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        --line_number;
                        break;
                    }

                    get_IV2D(trimmed_line);
                }
            }
            else if (key == "input_var_3d") {
                std::cout << "\nProcessing var_class:" << key << std::endl;
                while (true) {
                    std::streampos pos = file.tellg();
                    if (!std::getline(file, line)) {
                        break;
                    }
                    ++line_number;
                    trimmed_line = var_opf::trim(var_opf::get_key_before_char(line, '#'));

                    if (trimmed_line.empty()) {
                        continue;
                    }
                    key = var_opf::get_key_before_char(trimmed_line, ':');
                    if (key != trimmed_line) {
                        if (key.empty()) {
                            std::cerr << "Error: Illegal var_class name at line " << line_number << ": " << trimmed_line << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        if (!file.seekg(pos)) {
                            std::cerr << "Error: Failed to seek back in file" << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        --line_number;
                        break;
                    }

                    get_IV3D(trimmed_line);
                }
            }
            else if (key == "input_dimensionless_char_var") {
                std::cout << "\nProcessing var_class:" << key << std::endl;
                nondimensionalize_flag = true;
                while (true) {
                    std::streampos pos = file.tellg();
                    if (!std::getline(file, line)) {
                        break;
                    }
                    ++line_number;
                    trimmed_line = var_opf::trim(var_opf::get_key_before_char(line, '#'));

                    if (trimmed_line.empty()) {
                        continue;
                    }
                    key = var_opf::get_key_before_char(trimmed_line, ':');
                    if (key != trimmed_line) {
                        if (key.empty()) {
                            std::cerr << "Error: Illegal var_class name at line " << line_number << ": " << trimmed_line << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        if (!file.seekg(pos)) {
                            std::cerr << "Error: Failed to seek back in file" << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        --line_number;
                        break;
                    }

                    get_IDCV(trimmed_line);
                }
            }
            else if (key == "varying_bottom") {
                std::cout << "\nProcessing var_class:" << key << std::endl;
                varying_bottom_flag = true;
                while (true) {
                    std::streampos pos = file.tellg();
                    if (!std::getline(file, line)) {
                        break;
                    }
                    ++line_number;
                    trimmed_line = var_opf::trim(var_opf::get_key_before_char(line, '#'));

                    if (trimmed_line.empty()) {
                        continue;
                    }
                    key = var_opf::get_key_before_char(trimmed_line, ':');
                    if (key != trimmed_line) {

                        if (key.empty()) {
                            std::cerr << "Error: Illegal var_class name at line " << line_number << ": " << trimmed_line << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        // allow ':' exists in string muparser expression beta and gamma.
                        if (key.find('=') != std::string::npos) {
                            get_VB(trimmed_line);
                            continue;
                        }
                        else {
                            if (!file.seekg(pos)) {
                                std::cerr << "Error: Failed to seek back in file" << std::endl;
                                std::exit(EXIT_FAILURE);
                            }
                            --line_number;
                            break;
                        }
                        
                    }

                    get_VB(trimmed_line);
                }
            }
            else if (key == "save_data_control") {
                std::cout << "\nProcessing var_class:" << key << std::endl;
                while (true) {
                    std::streampos pos = file.tellg();
                    if (!std::getline(file, line)) {
                        break;
                    }
                    ++line_number;
                    trimmed_line = var_opf::trim(var_opf::get_key_before_char(line, '#'));

                    if (trimmed_line.empty()) {
                        continue;
                    }
                    key = var_opf::get_key_before_char(trimmed_line, ':');
                    if (key != trimmed_line) {
                        if (key.empty()) {
                            std::cerr << "Error: Illegal var_class name at line " << line_number << ": " << trimmed_line << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        if (!file.seekg(pos)) {
                            std::cerr << "Error: Failed to seek back in file" << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        --line_number;
                        break;
                    }

                    get_SDC(trimmed_line);
                }
            }
            else {
                std::cerr << "Error: Unknown var_class name: " << key << std::endl;
                std::exit(EXIT_FAILURE);
            }
            
        }
        std::cout << "\nEND get_IV()!\n" << std::endl;
    }

    bool check_IV() {
        std::cout << "\nStart  checking input variables !\n" << std::endl;
        if (GIV.dimension_num != 2 && GIV.dimension_num != 3) {
            std::cerr << "Error: dimension_num should be 2 or 3 in global_input_var " << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (GIV.DNO_M < 0 || GIV.DNO_Mb < 0) {
            std::cerr << "Error: DNO_M and DNO_Mb should be non-negative in global_input_var " << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (GIV.t_total == 0 || GIV.dt == 0) {
            std::cerr << "Error: t_total and dt should be non-zero in global_input_var " << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if ((GIV.t_total > 0 && GIV.dt < 0) || (GIV.t_total < 0 && GIV.dt > 0)) {
            std::cerr << "Error: t_total and dt should have the same sign in global_input_var " << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (std::fabs(GIV.t_total) < std::fabs(GIV.dt)) {
            std::cerr << "Error: abs(t_total) can't be smaller than abs(dt) in global_input_var " << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (GIV.mean_depth <= 0) {
            std::cerr << "Error: mean_depth should be positive in global_input_var " << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (GIV.gravity <= 0) {
            std::cerr << "Error: gravity should be positive in global_input_var " << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (GIV.NL_Ta<0) {
            std::cerr << "Error: NL_Ta should be non-negative in global_input_var " << std::endl;
            std::exit(EXIT_FAILURE);
        }

        // if (GIV.filter.empty()) presents not using filter.
        if (!GIV.filter.empty()) {
            bool is_exp = false;
            if (std::holds_alternative<std::string>(GIV.filter[0])) {
                if (std::get<std::string>(GIV.filter[0]) == "exp") {
                    is_exp = true;
                }
            }
            else {
                std::cerr << "Error: first item in global_input_var:filter should be a string(like exp)" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (is_exp) {
                if (GIV.filter.size() != 3) {
                    std::cerr << "Error: in global_input_var:filter, exp filter should have exactly two parameters after it (alpha, m)" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                
                auto is_positive_number = [](const auto& v) -> bool {
                    if (std::holds_alternative<int>(v)) {
                        return std::get<int>(v) > 0;
                    } else if (std::holds_alternative<double>(v)) {
                        return std::get<double>(v) > 0.0;
                    }
                    return false;
                };
                if (!is_positive_number(GIV.filter[1]) || !is_positive_number(GIV.filter[2])) {
                    std::cerr << "Error: in GIV.filter, exp filter parameters (alpha, m) must be positive numbers (int or double)" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            } // add other filter types check here...
            else {
                std::cerr << "Error: Unsupported filter type '";
                var_opf::print_variantIDBS_item(std::cerr, GIV.filter[0]);
                std::cerr << "' in GIV.filter." << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        
        if (GIV.dimension_num == 2) {
            if (IV2D.Lx <= 0) {
                std::cerr << "Error: Lx should be positive in input_var_2d " << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (IV2D.Nx <= 1) {
                std::cerr << "Error: Nx should be larger than 1 in input_var_2d " << std::endl;
                std::exit(EXIT_FAILURE);
            }
            
            if (!IV2D.initial_conditions.empty()) {
                bool is_Airy = false;
                if (std::holds_alternative<std::string>(IV2D.initial_conditions[0])) {
                    if (std::get<std::string>(IV2D.initial_conditions[0]) == "airy") {
                        is_Airy = true;
                    }// add more IC support here.
                }
                else {
                    std::cerr << "Error: first item in IV2D.initial_conditions should be a string(like Airy)" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                if (is_Airy) {
                    if (IV2D.initial_conditions.size() != 4) {
                        std::cerr << "Error: in IV2D.initial_conditions, Airy should have exactly three parameters after it (a(m), wave_length(m), theta_offset)" << std::endl;
                        std::exit(EXIT_FAILURE);
                    }
                    
                    auto is_positive_number = [](const auto& v) -> bool {
                        if (std::holds_alternative<int>(v)) {
                            return std::get<int>(v) > 0;
                        } else if (std::holds_alternative<double>(v)) {
                            return std::get<double>(v) > 0.0;
                        }
                        return false;
                    };
                    if (!is_positive_number(IV2D.initial_conditions[1]) || !is_positive_number(IV2D.initial_conditions[2])) {
                        std::cerr << "Error: in IV2D.initial_conditions, Airy a(m) and wave_length(m) parameters must be positive numbers (int or double)" << std::endl;
                        std::exit(EXIT_FAILURE);
                    }
                    else {
                        double Airy_a = var_opf::variantIDBS_to_double(IV2D.initial_conditions[1]);
                        if (Airy_a >= GIV.mean_depth) {
                            std::cerr << "Error: in IV2D.initial_conditions, Airy a(m) must be smaller than mean_depth." << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        else if (Airy_a < 0) {
                            std::cerr << "Error: in IV2D.initial_conditions, Airy a(m) must be positive." << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                        
                    }
                    if (!std::holds_alternative<int>(IV2D.initial_conditions[3]) && !std::holds_alternative<double>(IV2D.initial_conditions[3])) {
                        std::cerr << "Error: in IV2D.initial_conditions, Airy theta_offset parameter must be a number (int or double)" << std::endl;
                        std::exit(EXIT_FAILURE);
                    }
                } // add more IC check here.
                else {
                    std::cerr << "Error: Unsupported initial_conditions type '";
                    var_opf::print_variantIDBS_item(std::cerr, IV2D.initial_conditions[0]);
                    std::cerr << "' in IV2D.initial_conditions." << std::endl;
                    std::exit(EXIT_FAILURE);
                }

            }


            // check BC
            if (IV2D.boundary_conditions.empty()){
                std::cerr << "Error: IV2D.boundary_conditions can't be empty!  " << std::endl;
                std::exit(EXIT_FAILURE);
            }
            else {
                if (IV2D.boundary_conditions == "period") {
                    if (!IV2D.initial_conditions.empty() && std::get<std::string>(IV2D.initial_conditions[0]) == "airy") {
                        double n_of_L = IV2D.Lx / var_opf::variantIDBS_to_double(IV2D.initial_conditions[2]);
                        double distance = std::abs(n_of_L - std::round(n_of_L));
                        if (distance > 1e-10 ) {
                            std::cerr << "Error: IV2D.boundary_conditions if period , Lx(2D) should be an integer multiple of the Airy wavelength!  " << std::endl;
                            std::exit(EXIT_FAILURE);
                        }
                    }
                }
                else if (IV2D.boundary_conditions == "wall") {
                    mirror_extension_flag = true;
                }//add more BC check here.
                else {
                    std::cerr << "Error: IV2D.boundary_conditions must be period or wall!  " << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            }
        }
        else if (GIV.dimension_num == 3) { // to be implemented.
            if (IV3D.Lx <= 0 || IV3D.Ly <= 0) {
                std::cerr << "Error: Lx and Ly should be positive in input_var_3d !" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (IV3D.Nx <= 1 || IV3D.Ny <= 1) {
                std::cerr << "Error: Nx and Ny should be larger than 1 in input_var_3d !" << std::endl;
                std::exit(EXIT_FAILURE);
            }

            // check BC
            if (IV3D.boundary_conditions.empty()){
                std::cerr << "Error: IV3D.boundary_conditions can't be empty!  " << std::endl;
                std::exit(EXIT_FAILURE);
            }
            else {
                if (IV3D.boundary_conditions == "period"){
                    // add more check here.
                }
                else if (IV3D.boundary_conditions == "wall") {
                    mirror_extension_flag = true;
                }//add more BC check here.
                else {
                    std::cerr << "Error: IV3D.boundary_conditions must be period or wall!  " << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            }
        }
        
        if (nondimensionalize_flag){
            if (IDCV.L_char <= 0 || IDCV.T_char <= 0 || IDCV.M_char <= 0) {
                std::cerr << "Error: L_char T_char and M_char should be positive in input_dimensionless_char_var !" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }

        if (varying_bottom_flag) {
            // Ensure that VB.beta is a legal non-empty string and valid parser expression.
            if (VB.beta.empty()) {
                std::cerr << "VB.beta can't be empty!(should be an Parser expression) !" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (VB.time_dependency) {
                // bottom beta varying with time, should have \gamme = \partial_t \beta
                if (VB.gamma.empty()) {
                    std::cerr << "VB.gamma can't be empty if time_dependency!(should be an Parser expression) !" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            }
            mu::Parser parser_beta;
            mu::Parser parser_gamma;
            double x_temp = 0.0;
            double y_temp = 0.0;
            double t_temp = 0.0;
            if (GIV.dimension_num == 2) {
                parser_beta.DefineVar("x", &x_temp);
                parser_gamma.DefineVar("x", &x_temp);
                if (VB.time_dependency) {
                    parser_beta.DefineVar("t", &t_temp);
                    parser_gamma.DefineVar("t", &t_temp);
                }
            } else if (GIV.dimension_num == 3) {
                parser_beta.DefineVar("x", &x_temp);
                parser_beta.DefineVar("y", &y_temp);
                parser_gamma.DefineVar("x", &x_temp);
                parser_gamma.DefineVar("y", &y_temp);
                if (VB.time_dependency) {
                    parser_beta.DefineVar("t", &t_temp);
                    parser_gamma.DefineVar("t", &t_temp);
                }
            }

            try {
                parser_beta.SetExpr(VB.beta);
            } catch (mu::ParserError& e) {
                std::cerr << "VB.beta is an Illegal Parser expression !" << e.GetMsg() << std::endl;
                std::exit(EXIT_FAILURE);
            }

            if (VB.time_dependency) {
                try {
                    parser_gamma.SetExpr(VB.gamma);
                } catch (mu::ParserError& e) {
                    std::cerr << "VB.gamma is an Illegal Parser expression !" << e.GetMsg() << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            }
            

            mu::varmap_type beta_vars = parser_beta.GetUsedVar();

            if (!VB.time_dependency && beta_vars.find("t") != beta_vars.end()) {
                std::cerr << "Error: Expression contains 't' but time_dependency is false." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            else if (VB.time_dependency && beta_vars.find("t") == beta_vars.end()) {
                std::cerr << "Error: time_dependency is true but expression does not contain variable 't'." << std::endl;
                std::exit(EXIT_FAILURE);
            }

            if (GIV.dimension_num == 2 && beta_vars.find("y") != beta_vars.end()) {
                std::cerr << "Error: 2D simulation cannot use variable 'y'." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            
        }

        if (SDC.f_out <= 0) {
            std::cerr << "Error: in SDC.f_out, sample rate should be a positive number." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        if (SDC.out_file_path.empty()) {
            std::cerr << "Error: in SDC.out_file_path, out_file_path can't be empty!." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        var_opf::check_or_create_output_directory(SDC.out_file_path);

        std::cout << "\nEnd  checking input variables !\n" << std::endl;
        return true;
        
    }

    void nondimensionalize() {
        if (nondimensionalize_flag) {
            GIV.t_total /= IDCV.T_char;
            GIV.dt /= IDCV.T_char;
            GIV.mean_depth /= IDCV.L_char;
            GIV.gravity /= IDCV.L_char / (IDCV.T_char * IDCV.T_char);
            GIV.NL_Ta /= IDCV.T_char;
            if (GIV.dimension_num == 2) {
                IV2D.Lx /= IDCV.L_char;
                IV2D.x0 /= IDCV.L_char;
                if (!IV2D.initial_conditions.empty() && std::get<std::string>(IV2D.initial_conditions[0]) == "airy") {
                    var_opf::divide_at_indexofIDBS(IV2D.initial_conditions, IDCV.L_char , 1);
                    var_opf::divide_at_indexofIDBS(IV2D.initial_conditions, IDCV.L_char , 2);
                } // add other IC nondimensionalize here.
                
            }
            else if (GIV.dimension_num == 3) { // to be implemented!
                IV3D.Lx /= IDCV.L_char;
                IV3D.Ly /= IDCV.L_char;
                IV3D.x0 /= IDCV.L_char;
                IV3D.y0 /= IDCV.L_char;
            }

            // sampling rate nondimensionalize
            SDC.f_out *= IDCV.T_char;
            
        }
        else {
            std::cerr << "Error: Should not use nondimensionalize() if there's no var_class: input_dimensionless_char_var !" << std::endl;
            std::exit(EXIT_FAILURE); 
        }

    }

    void print_IV() {
        std::cout << std::fixed << std::setprecision(12);
        std::cout << "dimension_num = " << GIV.dimension_num << std::endl;
        std::cout << "DNO_M         = " << GIV.DNO_M << std::endl;
        std::cout << "DNO_Mb        = " << GIV.DNO_Mb << std::endl;
        std::cout << "Dealiase_flag = " << std::boolalpha << GIV.Dealiase_flag << std::endl;
        std::cout << "t_total       = " << GIV.t_total << std::endl;
        std::cout << "dt            = " << GIV.dt << std::endl;
        std::cout << "mean_depth    = " << GIV.mean_depth << std::endl;
        std::cout << "gravity       = " << GIV.gravity << std::endl;
        std::cout << "NL_Ta         = " << GIV.NL_Ta << std::endl;
        std::cout << "filter        = ";
        var_opf::print_variantIDBS_vector(std::cout, GIV.filter);
        std::cout << std::endl;

        std::cout << "2DLx          = " << IV2D.Lx << std::endl;
        std::cout << "2DNx          = " << IV2D.Nx << std::endl;
        std::cout << "2Dx0          = " << IV2D.x0 << std::endl;
        std::cout << "2D initial conditions        = ";
        var_opf::print_variantIDBS_vector(std::cout, IV2D.initial_conditions);
        std::cout << std::endl;
        std::cout << "2D boundary conditions      = " << IV2D.boundary_conditions << std::endl;


        std::cout << "3DLx          = " << IV3D.Lx << std::endl;
        std::cout << "3DLy          = " << IV3D.Ly << std::endl;
        std::cout << "3DNx          = " << IV3D.Nx << std::endl;
        std::cout << "3DNy          = " << IV3D.Ny << std::endl;
        std::cout << "3Dx0          = " << IV3D.x0 << std::endl;
        std::cout << "3Dy0          = " << IV3D.y0 << std::endl;
        std::cout << "3D initial conditions        = ";
        var_opf::print_variantIDBS_vector(std::cout, IV3D.initial_conditions);
        std::cout << std::endl;
        std::cout << "3D boundary conditions      = " << IV3D.boundary_conditions << std::endl;


        std::cout << "L_char        = " << IDCV.L_char << std::endl;
        std::cout << "T_char        = " << IDCV.T_char << std::endl;
        std::cout << "M_char        = " << IDCV.M_char << std::endl;

        
        std::cout << "beta          = " << VB.beta << std::endl;
        std::cout << "beta_time_dependency = " << VB.time_dependency << std::endl;
        std::cout << "gamma         = " << VB.gamma << std::endl;
        std::cout << "f_out         = " << SDC.f_out << std::endl;
        std::cout << "out_file_path = " << SDC.out_file_path << std::endl;


    }

} // namespace get_input_var
