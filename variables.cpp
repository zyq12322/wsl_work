#define _USE_MATH_DEFINES
#include "variables.h"
#include "var_op_function.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <stdexcept>




global_input_var GIV;
input_var_2d IV2D;
input_var_3d IV3D;
input_dimensionless_char_var IDCV;
varying_bottom VB;
save_data_control SDC;

bool nondimensionalize_flag = false;
bool varying_bottom_flag = false;
bool mirror_extension_flag = false;

internal_global_shared_var_2d IGSV2D;
internal_global_shared_var_3d IGSV3D;
internal_global_shared_temporary_var IGSTV;

std::unique_ptr<Dealiaser> dealiaser;


/**
 * @brief Writes a vector of variant values to an output file stream.
 *
 * Each element in the list is converted to its string representation using
 * `var_opf::variantIDBS_to_string()` and written to the file, separated by ", ".
 *
 * The output format is a single line of comma-separated values, e.g.:
 *     1, 2.5, true, hello
 *
 * @param ofs Output file stream to write to.
 * @param list Vector of variants (int, double, bool, std::string).
 */
void write_variantIDBS_list(
    std::ofstream& ofs,
    const std::vector<std::variant<int, double, bool, std::string>>& list)
{
    for (size_t i = 0; i < list.size(); ++i) {
        ofs << var_opf::variantIDBS_to_string(list[i]);
        if (i + 1 < list.size()) {
            ofs << ", ";
        }
    }
}

/**
 * @brief Writes the full simulation configuration to a text file.
 *
 * This function serializes multiple global configuration structures into a
 * human-readable configuration file format. The output is organized into
 * labeled sections depending on the simulation setup.
 *
 * The file includes:
 *
 * 1. Global simulation parameters (global_input_var)
 * 2. Dimension-dependent input variables (2D or 3D configuration)
 * 3. Optional nondimensionalization parameters
 * 4. Optional varying bottom parameters
 * 5. Output/save control settings
 *
 * The format is structured for readability and manual editing, using
 * key-value pairs and hierarchical section headers.
 *
 * @param filepath Path to the output configuration file.
 *
 * @throws std::runtime_error if the file cannot be opened.
 * @throws std::invalid_argument if dimension_num is not 2 or 3.
 */
void writeConfigFile(const std::string& filepath) {
    std::ofstream ofs(filepath);
    if (!ofs.is_open()) {
        throw std::runtime_error("Can't open file: " + filepath);
    }

    // set global output format
    ofs << std::defaultfloat << std::setprecision(15);

    // ========== 1. var_class : global_input_var ==========
    ofs << "global_input_var:\n";
    ofs << "  dimension_num = " << GIV.dimension_num << "\n";
    ofs << "  DNO_M = " << GIV.DNO_M << "\n";
    ofs << "  DNO_Mb = " << GIV.DNO_Mb << "\n";
    ofs << "  Dealiase_flag = " << std::boolalpha << GIV.Dealiase_flag << "\n";
    ofs << "  t_total = " << GIV.t_total << "\n";
    ofs << "  dt = " << GIV.dt << "\n";
    ofs << "  mean_depth = " << GIV.mean_depth << "\n";
    ofs << "  gravity = " << GIV.gravity << "\n";
    ofs << "  NL_Ta = " << GIV.NL_Ta << "\n";
    if (!GIV.filter.empty()){
        ofs << "  filter = ";
        write_variantIDBS_list(ofs, GIV.filter);
    }
    
    ofs << "\n\n";

    // ========== 2. select from IV2D and IV3D according to GIV.dimension_num ==========
    if (GIV.dimension_num == 2) {
        ofs << "input_var_2d:\n";
        ofs << "  Lx = " << IV2D.Lx << "\n";
        ofs << "  Nx = " << IV2D.Nx << "\n";
        ofs << "  x0 = " << IV2D.x0 << "\n";
        if (!IV2D.initial_conditions.empty()) {
            ofs << "  initial_conditions = ";
            write_variantIDBS_list(ofs, IV2D.initial_conditions);
            ofs << "\n";
        }
        
        ofs << "  boundary_conditions = " << IV2D.boundary_conditions << "\n\n";
    } else if (GIV.dimension_num == 3) {
        ofs << "input_var_3d:\n";
        ofs << "  Lx = " << IV3D.Lx << "\n";
        ofs << "  Ly = " << IV3D.Ly << "\n";
        ofs << "  Nx = " << IV3D.Nx << "\n";
        ofs << "  Ny = " << IV3D.Ny << "\n";
        ofs << "  x0 = " << IV3D.x0 << "\n";
        ofs << "  y0 = " << IV3D.y0 << "\n\n";
        if (!IV3D.initial_conditions.empty()) {
            ofs << "  initial_conditions = ";
            write_variantIDBS_list(ofs, IV3D.initial_conditions);
            ofs << "\n";
        }
        
        ofs << "  boundary_conditions = " << IV3D.boundary_conditions << "\n\n";
    } else {
        throw std::invalid_argument("dimension_num must be 2 or 3");
    }

    // ========== 3. optional：input_dimensionless_char_var ==========
    if (nondimensionalize_flag) {
        ofs << "input_dimensionless_char_var:\n";
        ofs << "  L_char = " << IDCV.L_char << "\n";
        ofs << "  T_char = " << IDCV.T_char << "\n";
        ofs << "  M_char = " << IDCV.M_char << "\n\n";
    }

    // ========== 4. optional：varying_bottom ==========
    if (varying_bottom_flag) {
        ofs << "varying_bottom:\n";
        ofs << "  time_dependency = " << std::boolalpha << VB.time_dependency << "\n";
        ofs << "  beta = " << VB.beta << "\n";
        ofs << "  gamma = " << VB.gamma << "\n\n";
    }

    // ========== 5. save_data_control ==========
    ofs << "save_data_control:\n";
    ofs << "  f_out = " << SDC.f_out << "\n";
    ofs << "  out_file_path = " << SDC.out_file_path << "\n";

    ofs.close();
}

// set IGSV2D, IGSV3D and IGSTV according to the IV(input vars)
void set_IGSV() {
    // choosing whether to perform mirror extension according to the boundary conditions
    // set IGSTV.N_space
    if (GIV.dimension_num == 2) {
        if (!mirror_extension_flag) { // periodic BC
            IGSTV.N_space = static_cast<size_t>(IV2D.Nx);
            
        }
        else { // mirror extension
            IGSTV.N_space = static_cast<size_t>(IV2D.Nx * 2);
        }
    }
    else if (GIV.dimension_num == 3) {
        if (!mirror_extension_flag) { // periodic BC
            IGSTV.N_space = static_cast<size_t>(IV3D.Nx * IV3D.Ny);
        }
        else { // mirror extension
            IGSTV.N_space = static_cast<size_t>(IV3D.Nx * IV3D.Ny * 4);
        }
    }
    // set IGSTV.t_count
    // if start from t=0;
    IGSTV.t_count = 0.0;

    // set IGSTV.sampling_interval
    IGSTV.sampling_interval = 1/ SDC.f_out;

    // if var_class varying_bottom is inputed.
    // then set IGSTV.beta_parser and gamma_parser
    if (varying_bottom_flag) {
        IGSTV.beta_parser.SetExpr(VB.beta);
        if (VB.time_dependency) {
            IGSTV.gamma_parser.SetExpr(VB.gamma);
        }
        if (GIV.dimension_num == 2) {
            IGSTV.beta_parser.DefineVar("x", &IGSTV.vb_var_x);
            IGSTV.gamma_parser.DefineVar("x", &IGSTV.vb_var_x);
            if (VB.time_dependency) {
               IGSTV.beta_parser.DefineVar("t", &IGSTV.vb_var_t); 
               IGSTV.gamma_parser.DefineVar("t", &IGSTV.vb_var_t);
            }
        }
        else if (GIV.dimension_num == 3) {
            IGSTV.beta_parser.DefineVar("x", &IGSTV.vb_var_x);
            IGSTV.beta_parser.DefineVar("y", &IGSTV.vb_var_y);
            IGSTV.gamma_parser.DefineVar("x", &IGSTV.vb_var_x);
            IGSTV.gamma_parser.DefineVar("y", &IGSTV.vb_var_y);
            if (VB.time_dependency) {
               IGSTV.beta_parser.DefineVar("t", &IGSTV.vb_var_t); 
               IGSTV.gamma_parser.DefineVar("t", &IGSTV.vb_var_t); 
            }
        }
    }
    
    // alloc memory for vectors of free-surface quantities in IGSTV.
    resize_IGSTV_fs_quantities(GIV.DNO_M, GIV.DNO_Mb, IGSTV.N_space);

    // set IGSV2D
    if (GIV.dimension_num == 2) {
        // set IGSV2D.x_axis
        if (!mirror_extension_flag) {
            // dimensionless if nondimensionalize_flag
            IGSV2D.x_axis = var_opf::linspace(IV2D.x0, IV2D.x0 + IV2D.Lx, IGSTV.N_space, false);
        }
        else {
            //mirror extension part [0, IV2D.Nx-1]
            // dimensionless if nondimensionalize_flag
            IGSV2D.x_axis = var_opf::linspace(IV2D.x0 - IV2D.Lx , IV2D.x0 + IV2D.Lx, IGSTV.N_space, false);
        }
        // set IGSV2D.dx
        IGSV2D.dx = IV2D.Lx / IV2D.Nx;
        // set IGSV2D.k_axis
        IGSV2D.k_axis = var_opf::fftfreq(IGSTV.N_space, IGSV2D.dx);
        // 1/length - > k(wave_number)
        for (double& k : IGSV2D.k_axis) {
            k *= 2.0 * M_PI;
            //std::cout << k << std::endl;
        }

        // alloc memory for vectors in IGSV2D
        IGSV2D.tanh_hD.resize(IGSTV.N_space);
        IGSV2D.sech_hD.resize(IGSTV.N_space);
        IGSV2D.G00.resize(IGSTV.N_space);

        int max_order = std::max({GIV.DNO_M, GIV.DNO_Mb, 1});
        IGSV2D.D_pow_list.resize(max_order + 1);
        for (int p = 0; p <= max_order; ++p) {
            IGSV2D.D_pow_list[p].resize(IGSTV.N_space);
        }
        IGSV2D.fact_list.resize(max_order + 1);
        
        // set vectors in IGSV2D
        std::transform(IGSV2D.k_axis.begin(), IGSV2D.k_axis.end(),
                    IGSV2D.tanh_hD.begin(),
                    [h = GIV.mean_depth](double k) { return std::tanh(k * h); });
            
        std::transform(IGSV2D.k_axis.begin(), IGSV2D.k_axis.end(),
                    IGSV2D.sech_hD.begin(),
                    [h = GIV.mean_depth](double k) { return 1.0 / std::cosh(k * h); });

            
        for (size_t i = 0; i < IGSTV.N_space; ++i) {
            IGSV2D.G00[i] = IGSV2D.tanh_hD[i] * IGSV2D.k_axis[i];
        }

        std::fill(IGSV2D.D_pow_list[0].begin(), IGSV2D.D_pow_list[0].end(), 1.0);

        IGSV2D.D_pow_list[1] = IGSV2D.k_axis;

        for (int p = 2; p <= max_order; ++p) {
            const auto& prev = IGSV2D.D_pow_list[p-1];
            auto& curr = IGSV2D.D_pow_list[p];
            for (size_t i = 0; i < IGSTV.N_space; ++i) {
                curr[i] = prev[i] * IGSV2D.k_axis[i];
            }
        }
        
        IGSV2D.fact_list[0] = 1.0;
        for (int i = 1; i <= max_order; ++i) {
            IGSV2D.fact_list[i] = IGSV2D.fact_list[i-1] * static_cast<double>(i);
            //std::cout << IGSV2D.fact_list[i] << std::endl;
        }

    }
}
/**
 * @brief Sets the 1D varying bottom (beta) profile and its FFT for a 2D problem.
 *
 * This function computes the spatial distribution of the bottom topography
 * (or any x‑dependent coefficient) at time `t`(dimensionless if nondimensionalize_flag). 
 * The result is stored in
 * `IGSTV.beta` (real vector) and its forward FFT in `IGSTV.beta_hat`.
 *
 * The calculation respects:
 * - Dimensionality: if `nondimensionalize_flag` is true, all inputs to the
 *   parser are multiplied by characteristic scales and the output is divided
 *   by the length scale; otherwise raw values are used.
 * - Boundary conditions:
 *   - "period": beta is defined on [0, Nx-1] and used as‑is.
 *   - "wall": the domain is mirrored to create a symmetric extension of length
 *     2*Nx, suitable for sine/cosine transforms or wall‑bounded problems.
 *
 * @param t   Current time (dimensional or dimensionless according to flag).
 *
 * @pre   `varying_bottom_flag` must be true.
 * @pre   `IGSTV.beta_parser` is properly configured.
 * @pre   `IGSV2D.x_axis` and `IV2D.Nx` are correctly set.
 *
 * @post  `IGSTV.beta` contains the real profile.
 * @post  `IGSTV.beta_hat` contains its 1D complex FFT (size matches beta).
 *
 * @note  The FFT plan is fixed, make sure IGSTV.N_space not change.
 *
 * @warning The function exits the program if `varying_bottom_flag` is false.
 *
 * @bug    Mirror extension part 0,1,..., IV2D.Nx-1
 */
// input and output dimensionless if nondimensionalize_flag, else dimensional.
void set_beta_2d(double t) {// set IGSTV.beta and IGSTV.beta_hat
    static var_opf::fft1d fft(IGSTV.N_space);
    static double L_scale = nondimensionalize_flag ? IDCV.L_char: 1.0;
    static double T_scale = nondimensionalize_flag ? IDCV.T_char: 1.0;
    // set time IGSTV.vb_var_t with dimension.
    if (varying_bottom_flag) {
        IGSTV.vb_var_t = t * T_scale;

        if (!mirror_extension_flag) {
            for (int i =0; i< IV2D.Nx; ++i) {
                IGSTV.vb_var_x = IGSV2D.x_axis[i] * L_scale;
                IGSTV.beta[i] = IGSTV.beta_parser.Eval() / L_scale;
            }
        }
        else {//mirror extension
            for (int i = IV2D.Nx; i< IV2D.Nx * 2; ++i) {
                IGSTV.vb_var_x = IGSV2D.x_axis[i] * L_scale;
                IGSTV.beta[i] = IGSTV.beta_parser.Eval() / L_scale;
                IGSTV.beta[IV2D.Nx * 2 - i] = IGSTV.beta[i];
            }
            // linear extension
            IGSTV.beta[0] = 2* IGSTV.beta[1] - IGSTV.beta[2];
        }
        
        IGSTV.beta_hat = fft.fft(var_opf::double_to_complex_real(IGSTV.beta));
    }
    else {
        std::cerr << "Error: use set_beta() without inputing varclass:varying_bottom!"<< std::endl;
        std::exit(EXIT_FAILURE);
    }
}

/**
 * @brief Evaluates the 1D gamma profile (time‑varying bottom‑related field) and its FFT for a 2D problem.
 *
 * This function computes the spatial distribution of gamma (e.g., a bottom friction coefficient
 * or sediment transport rate) at time `t`. The result is stored in `IGSTV.gamma` (real vector)
 * and its forward FFT in `IGSTV.gamma_hat`.
 *
 * The calculation respects:
 * - Dimensionality: if `nondimensionalize_flag` is true, time `t` is multiplied by `IDCV.T_char`,
 *   x‑coordinates are multiplied by `IDCV.L_char`, and the parser output is divided by `IDCV.L_char`
 *   and multiplied by `IDCV.T_char` (to obtain correct dimensions for gamma).
 * - Boundary conditions:
 *   - "period": gamma is defined on indices [0, Nx-1] (length Nx).
 *   - "wall": symmetric mirror extension to length 2*Nx. Indices [Nx, 2Nx-1] are evaluated,
 *             then mirrored to indices [1, Nx-1]; index 0 is linearly extrapolated.
 *
 * @param t   Current time. Dimensionless if `nondimensionalize_flag` is true, else dimensional.
 *
 * @pre   `varying_bottom_flag` and `VB.time_dependency` must be true.
 * @pre   `IGSTV.gamma_parser` is properly configured and uses `IGSTV.gamma_var_t`/`gamma_var_x`
 *        (or a shared `t_var`/`x_var`) for evaluation.
 * @pre   `IV2D.Nx`, `IGSV2D.x_axis` are correctly set. For "wall" BC, `x_axis` and `IGSTV.gamma`
 *        must have size at least 2*Nx.
 *
 * @post  `IGSTV.gamma` contains the profile (size Nx for "period", 2*Nx for "wall").
 * @post  `IGSTV.gamma_hat` contains its 1D complex FFT (same size).
 *
 * @note  The FFT plan is created once using `IGSTV.N_space` and reused. `IGSTV.N_space`
 *        is assumed constant after the first call.
 * @note  The scaling factors `L_scale` and `T_scale` are evaluated only at first call;
 *        changing `nondimensionalize_flag` later will have no effect.
 *
 * @warning The function exits the program if the boundary condition is unknown or
 *          the varying‑bottom flags are false.
 *
 * @bug    The function currently writes to `IGSTV.beta_var_t` and `IGSTV.beta_var_x` instead
 *        of dedicated gamma variables. This can corrupt beta state and should be fixed.
 */
// input and output dimensionless if nondimensionalize_flag, else dimensional.
void set_gamma_2d(double t) {// set IGSTV.gamma and IGSTV.gamma_hat
    static var_opf::fft1d fft(IGSTV.N_space);
    static double L_scale = nondimensionalize_flag ? IDCV.L_char: 1.0;
    static double T_scale = nondimensionalize_flag ? IDCV.T_char: 1.0;
    if (varying_bottom_flag && VB.time_dependency) {
        
        IGSTV.vb_var_t = t * T_scale;

        if (!mirror_extension_flag) {
            for (int i =0; i< IV2D.Nx; ++i) {
                IGSTV.vb_var_x = IGSV2D.x_axis[i] * L_scale;
                IGSTV.gamma[i] = IGSTV.gamma_parser.Eval() / L_scale * T_scale;
            }
        }
        else {//mirror extension
            for (int i = (IV2D.Nx); i< IV2D.Nx * 2; ++i) {
                IGSTV.vb_var_x = IGSV2D.x_axis[i] * L_scale;
                IGSTV.gamma[i] = IGSTV.gamma_parser.Eval() / L_scale * T_scale;
                IGSTV.gamma[IV2D.Nx * 2 - i] = IGSTV.gamma[i];
            }
            IGSTV.gamma[0] = 2* IGSTV.gamma[1] - IGSTV.gamma[2];
        }
            
        IGSTV.gamma_hat = fft.fft(var_opf::double_to_complex_real(IGSTV.gamma));
    }
    else {
        std::cerr << "Error: use set_gamma() without inputing varclass:varying_bottom or bottom doesn't vary with time!"<< std::endl;
        std::exit(EXIT_FAILURE);
    }
}
