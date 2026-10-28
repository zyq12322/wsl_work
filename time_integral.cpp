#define _USE_MATH_DEFINES
#include "time_integral.h"
#include "var_op_function.h"
#include "variables.h"
#include <cmath>
#include <algorithm>
#include <cstdlib>

namespace time_integral {

    void initial_conditions() {// please call this after nondimensionalize
        
        if (GIV.dimension_num == 2) {

            // check the size of 2d arrays.
            if (IGSTV.eta.size() != IGSV2D.x_axis.size()) {
                std::cerr << "Error: eta size (" << IGSTV.eta.size() 
                        << ") is not equal to x_axis size (" << IGSV2D.x_axis.size() 
                        << ") in initial_conditions." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (IGSTV.phis.size() != IGSV2D.x_axis.size()) {
                std::cerr << "Error: phis size (" << IGSTV.phis.size() 
                        << ") is not equal to x_axis size (" << IGSV2D.x_axis.size() 
                        << ") in initial_conditions." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (IGSTV.beta.size() != IGSV2D.x_axis.size()) {
                std::cerr << "Error: beta size (" << IGSTV.beta.size() 
                        << ") is not equal to x_axis size (" << IGSV2D.x_axis.size() 
                        << ") in initial_conditions." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (IGSTV.gamma.size() != IGSV2D.x_axis.size()) {
                std::cerr << "Error: gamma size (" << IGSTV.gamma.size() 
                        << ") is not equal to x_axis size (" << IGSV2D.x_axis.size() 
                        << ") in initial_conditions." << std::endl;
                std::exit(EXIT_FAILURE);
            }


            // initial_conditions for eta and phis-in both physical space and Fourier space;
            if (IV2D.initial_conditions.empty()) {// all zeros IC
                std::fill(IGSTV.eta.begin(), IGSTV.eta.end(), 0.0);
                std::fill(IGSTV.phis.begin(), IGSTV.phis.end(), 0.0);
            }
            else if (std::get<std::string>(IV2D.initial_conditions[0]) == "airy") {
                double a_val = var_opf::variantIDBS_to_double(IV2D.initial_conditions[1]);
                double l_val = var_opf::variantIDBS_to_double(IV2D.initial_conditions[2]);
                double off_val = var_opf::variantIDBS_to_double(IV2D.initial_conditions[3]);
                double k_val = 2 * M_PI / l_val;
                double h = GIV.mean_depth;
                double omega = std::sqrt(GIV.gravity * k_val * std::tanh(k_val * h) );
                if (!mirror_extension_flag) {
                    for (int i = 0; i < IV2D.Nx; ++i) {
                        IGSTV.eta[i] = a_val * std::cos( k_val * IGSV2D.x_axis[i] + off_val);
                        IGSTV.phis[i] = omega * a_val * std::cosh(k_val * (h + IGSTV.eta[i]))
                        / k_val / std::sinh(k_val * h) * std::sin(k_val * IGSV2D.x_axis[i] + off_val);
                    }
                }
                else {
                    for (int i = IV2D.Nx; i < IV2D.Nx * 2; ++i) {
                        IGSTV.eta[i] = a_val * std::cos( k_val * IGSV2D.x_axis[i] + off_val);
                        IGSTV.eta[IV2D.Nx * 2 - i] = IGSTV.eta[i];
                    }
                    IGSTV.eta[0] = 2* IGSTV.eta[1] - IGSTV.eta[2];
                    std::fill(IGSTV.phis.begin(), IGSTV.phis.end(), 0.0);
                }
                
            }
            else {
                std::cerr << "Error: Unsupported initial_conditions type '";
                var_opf::print_variantIDBS_item(std::cerr, IV2D.initial_conditions[0]);
                std::cerr << "' in IV2D.initial_conditions." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            IGSTV.eta_hat = fft(IGSTV.eta);
            IGSTV.phis_hat = fft(IGSTV.phis);

            // initial_conditions for beta and beta_hat gamma, gamma_hat -in both physical space and Fourier space;
            if (varying_bottom_flag) {
                // set IGSTV.beta and IGSTV.beta_hat
                set_beta_2d(IGSTV.t_count);

                // set IGSTV.gamma and IGSTV.gamma_hat
                if (VB.time_dependency) {
                    set_gamma_2d(IGSTV.t_count);
                }
                else {
                    std::complex<double> zero_complex(0.0,0.0);
                    std::fill(IGSTV.gamma.begin(), IGSTV.gamma.end(), 0.0);
                    std::fill(IGSTV.gamma_hat.begin(), IGSTV.gamma_hat.end(), zero_complex);
                }
            }
            else {
                std::complex<double> zero_complex(0.0,0.0);

                // set IGSTV.beta and IGSTV.beta_hat
                std::fill(IGSTV.beta.begin(), IGSTV.beta.end(), 0.0);
                std::fill(IGSTV.beta_hat.begin(), IGSTV.beta_hat.end(), zero_complex);

                // set IGSTV.gamma and IGSTV.gamma_hat
                std::fill(IGSTV.gamma.begin(), IGSTV.gamma.end(), 0.0);
                std::fill(IGSTV.gamma_hat.begin(), IGSTV.gamma_hat.end(), zero_complex);
            }

        }
        else if (GIV.dimension_num == 3) {
            std::cerr << "Error: 3D simulation not yet implemented.\n";
            std::exit(EXIT_FAILURE);
        }
        else {
            std::cerr << "Error: dimension_num should be 2 or 3 in global_input_var " << std::endl;
            std::exit(EXIT_FAILURE);
        }

    }

}

