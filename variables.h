#ifndef VAR_H
#define VAR_H

#include "dealiaser.h"
#include <vector>
#include <string>
#include <variant>
#include <complex>
#include <muParser.h>
#include <fstream>



//input vars

// GIV
struct global_input_var{
    int dimension_num = 2;
    int DNO_M = 4;
    int DNO_Mb = 6;
    bool Dealiase_flag = true;
    double t_total = 0.0;//dimensionless if nondimensionalize_flag;
    double dt = 0.0;//dimensionless if nondimensionalize_flag;
    double mean_depth = 0.0;//dimensionless if nondimensionalize_flag;
    double gravity = 9.81;//dimensionless if nondimensionalize_flag;
    std::vector<std::variant<int, double, bool, std::string>> filter;
    double NL_Ta = 0.0;// dimensionless if nondimensionalize_flag
};
// IV2D
struct input_var_2d{
    double Lx = 0;//dimensionless if nondimensionalize_flag;
    int Nx = 0;
    double x0 = 0;//dimensionless if nondimensionalize_flag;
    std::vector<std::variant<int, double, bool, std::string>> initial_conditions;//dimensionless if nondimensionalize_flag;
    std::string boundary_conditions;
};
// IV3D
struct input_var_3d{
    double Lx = 0;//dimensionless if nondimensionalize_flag;
    double Ly = 0;//dimensionless if nondimensionalize_flag;
    int Nx = 0;
    int Ny = 0;
    double x0 = 0;//dimensionless if nondimensionalize_flag;
    double y0 = 0;//dimensionless if nondimensionalize_flag;
    std::vector<std::variant<int, double, bool, std::string>> initial_conditions;//dimensionless if nondimensionalize_flag;
    std::string boundary_conditions;
};
// IDCV
struct input_dimensionless_char_var{
    double L_char = 0;// m
    double T_char = 0;// s
    double M_char = 0;// kg
};
// VB
struct varying_bottom{
    // if (!time_dependency) ignore string gamma
    bool time_dependency = false;
    // with dimension-x(m),y(m),t(s),beta(m)
    std::string beta;
    // with dimension-x(m),y(m),t(s),gamma(m/s)
    std::string gamma; // \partial_t \beta
};
// SDC
struct save_data_control {
    double f_out = 0;//dimensionless if nondimensionalize_flag;
    std::string out_file_path;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// internal vars

// IGSV2D<-set_IGSV()
struct internal_global_shared_var_2d {
    //all params are dimensionless if nondimensionalize_flag;
    std::vector<double> x_axis;//len == Nx(BC:period),2*Nx(BC:wall);
    double dx;//dimensionless if nondimensionalize_flag;
    
    std::vector<double> k_axis;//len == Nx(BC:period),2*Nx(BC:wall);
    std::vector<double> tanh_hD;//len == Nx(BC:period),2*Nx(BC:wall);
    std::vector<double> sech_hD;//len == Nx(BC:period),2*Nx(BC:wall);
    std::vector<double> G00;//len == Nx(BC:period),2*Nx(BC:wall);
    std::vector<std::vector<double>> D_pow_list;//len == std::max({GIV.DNO_M, GIV.DNO_Mb, 1}) + 1
    std::vector<double> fact_list;
};
// IGSV3D<-set_IGSV()
struct internal_global_shared_var_3d {
    //all params are dimensionless if nondimensionalize_flag;
    std::vector<double> x_axis;// len == Nx(BC:period),2*Nx(BC:wall);
    std::vector<double> y_axis;// len == Ny(BC:period),2*Ny(BC:wall);
    double dx, dy;
    
    std::vector<double> kx_axis;// len == Nx(BC:period),2*Nx(BC:wall);
    std::vector<double> ky_axis;// len == Ny(BC:period),2*Ny(BC:wall);
    std::vector<double> tanh_h_absD;// len == Nx*Ny(BC:period),4*Nx*Ny(BC:wall);
    std::vector<double> sech_h_absD;// len == Nx*Ny(BC:period),4*Nx*Ny(BC:wall);
    std::vector<double> G00; // len == Nx*Ny(BC:period),4*Nx*Ny(BC:wall);
    std::vector<std::vector<double>> D_pow_list;
    std::vector<double> fact_list;
};

// IGSTV<-set_IGSV()
// shared global variables-if nondimensionalize_flag:dimensionless;else: dimensional
struct internal_global_shared_temporary_var {
    // total points number in space.
    size_t N_space = 0;// 2d: == Nx(BC:period),2*Nx(BC:wall); 3d: == Nx*Ny(BC:period),4*Nx*Ny(BC:wall);
    double t_count = 0.0;//if nondimensionalize_flag:dimensionless;else: dimensional
    // store sampling time interval-if nondimensionalize_flag:dimensionless;else: dimensional
    double sampling_interval = 0.0;// if nondimensionalize_flag:dimensionless;
    // parser expression for varying bottom.(with dimension)
    mu::Parser beta_parser;
    mu::Parser gamma_parser;
    // varying bottom arguments for parser string expression beta_parser and gamma_parser binding.
    // with dimension
    double vb_var_t = 0.0;
    double vb_var_x = 0.0;
    double vb_var_y = 0.0;
    
    // if nondimensionalize_flag:dimensionless;
    // physical quantities on the free surface
    std::vector<double> eta;// 2d:len == Nx(BC:period),2*Nx(BC:wall); 3d:len == Nx*Ny(BC:period),4*Nx*Ny(BC:wall);
    std::vector<double> phis;// 2d:len == Nx(BC:period),2*Nx(BC:wall); 3d:len == Nx*Ny(BC:period),4*Nx*Ny(BC:wall);
    std::vector<double> beta;// 2d:len == Nx(BC:period),2*Nx(BC:wall); 3d:len == Nx*Ny(BC:period),4*Nx*Ny(BC:wall);
    std::vector<double> gamma;// 2d:len == Nx(BC:period),2*Nx(BC:wall); 3d:len == Nx*Ny(BC:period),4*Nx*Ny(BC:wall);
    // physical quantities on the free surface in Fourier space (wavenumber axis)
    std::vector<std::complex<double>> eta_hat;// 2d:len == Nx(period),2*Nx(wall); 3d:len == Nx*Ny(BC:period),4*Nx*Ny(BC:wall);
    std::vector<std::complex<double>> phis_hat;// 2d:len == Nx(period),2*Nx(wall); 3d:len == Nx*Ny(BC:period),4*Nx*Ny(BC:wall);
    std::vector<std::complex<double>> beta_hat;// 2d:len == Nx(period),2*Nx(wall); 3d:len == Nx*Ny(BC:period),4*Nx*Ny(BC:wall);
    std::vector<std::complex<double>> gamma_hat;// 2d:len == Nx(period),2*Nx(wall); 3d:len == Nx*Ny(BC:period),4*Nx*Ny(BC:wall);
    // pow list in Fourier space.-for DNO_calculation.
    std::vector<std::vector<std::complex<double>>> eta_pow_hat;// (DNO_M + 1, eta.size())
    std::vector<std::vector<std::complex<double>>> beta_pow_hat;// (DNO_Mb + 1, beta.size())
};


// declare global variables.
// store input var
extern global_input_var GIV;
extern input_var_2d IV2D;
extern input_var_3d IV3D;
extern input_dimensionless_char_var IDCV;
extern varying_bottom VB;
extern save_data_control SDC;

// section exits flags
// an indicator of whether var_class input_dimensionless_char_var is present.
extern bool nondimensionalize_flag; // set in get_input_var::get_IV(input_filename);
// an indicator of whether var_class varying_bottom is present.
extern bool varying_bottom_flag;// set in get_input_var::get_IV(input_filename);

// mirror extension flag.
// use mirror extension or not.
extern bool mirror_extension_flag;

// internal shared variables
extern internal_global_shared_var_2d IGSV2D;
extern internal_global_shared_var_3d IGSV3D;
extern internal_global_shared_temporary_var IGSTV;

// internal shared tools.
// dealiaser for dealiasing product-global class case.
extern std::unique_ptr<Dealiaser> dealiaser;


void write_variantIDBS_list(
    std::ofstream& ofs,
    const std::vector<std::variant<int, double, bool, std::string>>& list);

void writeConfigFile(const std::string& filepath);

/**
 * @brief Resizes all primary and spectral field storage arrays.
 *
 * This function prepares the internal simulation storage (IGSTV) by resizing
 * both physical-space fields and spectral-space (Fourier) fields to match
 * the current spatial resolution, as well as allocating power-series spectral
 * buffers.
 *
 * Specifically, it performs the following operations:
 *
 * 1. Resizes physical-space fields to size N_space:
 *    - eta, phis, beta, gamma
 *
 * 2. Resizes corresponding spectral (hat) fields:
 *    - eta_hat, phis_hat, beta_hat, gamma_hat
 *
 * 3. Allocates power-series spectral storage:
 *    - eta_pow_hat with (DNO_M + 1) levels
 *    - beta_pow_hat with (DNO_Mb + 1) levels
 *    - Each level is resized to N_space
 *
 * This ensures all required arrays are properly initialized before any
 * numerical computation or FFT operations are performed.
 *
 * @param DNO_M Maximum polynomial/spectral order for eta powers.
 * @param DNO_Mb Maximum polynomial/spectral order for beta powers.
 * @param N_space Number of spatial grid points.
 */
inline void resize_IGSTV_fs_quantities(int DNO_M , int DNO_Mb , size_t N_space) {
    IGSTV.eta.resize(N_space);
    IGSTV.phis.resize(N_space);
    IGSTV.beta.resize(N_space);
    IGSTV.gamma.resize(N_space);
    IGSTV.eta_hat.resize(N_space);
    IGSTV.phis_hat.resize(N_space);
    IGSTV.beta_hat.resize(N_space);
    IGSTV.gamma_hat.resize(N_space);
    IGSTV.eta_pow_hat.resize(DNO_M + 1);
    IGSTV.beta_pow_hat.resize(DNO_Mb + 1);
    for (int p = 0; p <= DNO_M; ++p) {
        IGSTV.eta_pow_hat[p].resize(N_space);
    }
    for (int p = 0; p <= DNO_Mb; ++p) {
        IGSTV.beta_pow_hat[p].resize(N_space);
    }
}

// set IGSV2D or IGSV3D(according to GIV.dimension_num) , 
// set scalars in IGSTV.
// alloc memory for vectors of all free-surface quantities in IGSTV.
// use dimensionalized variables if nondimensionalize_flag.
void set_IGSV();
// set IGSTV.beta and IGSTV.beta_hat
// use if GIV.dimension_num == 2 and varying_bottom_flag
// input dimensionless t if nondimensionalize_flag
// set dimensionless beta if nondimensionalize_flag
void set_beta_2d(double t);
// set IGSTV.gamma and IGSTV.gamma_hat
// use if GIV.dimension_num == 2 and varying_bottom_flag && VB.time_dependency
// input dimensionless t if nondimensionalize_flag
// set dimensionless gamma if nondimensionalize_flag
void set_gamma_2d(double t);
#endif