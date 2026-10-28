#ifndef GETIV_H
#define GETIV_H


#include <string>

namespace process_input_var {

    void get_GIV(const std::string& trimmed_line);
    void get_IV2D(const std::string& trimmed_line);
    void get_IV3D(const std::string& trimmed_line);
    void get_IDCV(const std::string& trimmed_line);
    void get_VB(const std::string& trimmed_line);
    void get_SDC(const std::string& trimmed_line);


    void get_IV(const std::string& input_file);
    bool check_IV();
    void nondimensionalize();
    void print_IV();
} // namespace get_input_var



#endif