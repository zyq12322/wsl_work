#ifndef DNO_CALC_H
#define DNO_CALC_H

#include "dealiaser.h"
#include "variables.h"
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <cassert>

void G1_m0_0_2d(
    const std::vector<std::complex<double>>& x_hat,//len == Nx(period), 2*Nx(wall)
    const std::vector<std::vector<std::complex<double>>>& G2_m0_m1_list,
    std::vector<std::vector<std::complex<double>>>& G1_m0_m1_list,
    int DNO_M = GIV.DNO_M,
    int DNO_Mb = GIV.DNO_Mb);//fill this with G_1^{(m0,0)},m0=0,1,...,DNO_Mb
void G2_m0_0_2d(
    const std::vector<std::complex<double>>& x_hat,//len == Nx(period), 2*Nx(wall)
    std::vector<std::vector<std::complex<double>>>& G2_m0_m1_list,
    int DNO_M = GIV.DNO_M,
    int DNO_Mb = GIV.DNO_Mb);//fill this with G_2^{(m0,0)},m0=0,1,...,DNO_Mb

void G1_0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
    const std::vector<std::complex<double>>& x_hat,//len == Nx(period), 2*Nx(wall)
    std::vector<std::vector<std::complex<double>>>& G1_m0_m1_list,
    int DNO_M = GIV.DNO_M);//fill this with G_1^{(0,m1)},m1=0,1,...,DNO_M

void G2_0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
    const std::vector<std::complex<double>>& x_hat,//len == Nx(period), 2*Nx(wall)
    std::vector<std::vector<std::complex<double>>>& G1_m0_m1_list,
    std::vector<std::vector<std::complex<double>>>& G2_m0_m1_list,//fill this with G_2^{(0,m1)},m1=0,1,...,DNO_M
    int DNO_M = GIV.DNO_M);

void G1_m0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
    int m0,
    const std::vector<std::complex<double>>& x_hat,//len == Nx(period), 2*Nx(wall)
    std::vector<std::vector<std::complex<double>>>& G2_m0_m1_list,
    std::vector<std::vector<std::complex<double>>>& G1_m0_m1_list,
    int DNO_M = GIV.DNO_M,
    int DNO_Mb = GIV.DNO_Mb);//fill this with G_1^{(m0,m1)},m1=0,1,...,DNO_M

void G2_m0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
    int m0,
    const std::vector<std::complex<double>>& x_hat,//len == Nx(period), 2*Nx(wall)
    std::vector<std::vector<std::complex<double>>>& G1_m0_m1_list,
    std::vector<std::vector<std::complex<double>>>& G2_m0_m1_list,
    int DNO_M = GIV.DNO_M,
    int DNO_Mb = GIV.DNO_Mb);//fill this with G_2^{(m0,m1)},m1=1,...,DNO_M

inline void resize_G_m0_m1_list_2d(std::vector<std::vector<std::complex<double>>>& Gi_m0_m1_list, 
                                    int DNO_Mb,
                                    int DNO_M,
                                    size_t x_size) {
    Gi_m0_m1_list.resize((DNO_Mb+1) * (DNO_M+1));
    for (int i =0;i < (DNO_Mb+1) * (DNO_M+1); ++i) {
        Gi_m0_m1_list[i].resize(x_size);
    }
}
//
inline void set_eta_beta_pow_hat(Dealiaser1D& dealiaser, 
                                const std::vector<std::complex<double>>& eta_hat, 
                                const std::vector<std::complex<double>>& beta_hat,
                                int DNO_Mb = GIV.DNO_Mb,
                                int DNO_M = GIV.DNO_M)
{
    std::fill(IGSTV.eta_pow_hat[0].begin(), IGSTV.eta_pow_hat[0].end(), std::complex<double>(1.0, 0.0));
    if (DNO_M > 0) {
        IGSTV.eta_pow_hat[1] = eta_hat;
    }
    
    for (int p = 2; p <= DNO_M; ++p) {
        if (GIV.Dealiase_flag) {
            dealiaser(IGSTV.eta_pow_hat[p-1].data(), eta_hat.data(),
                    false,   // input in Fourier space
                    false,   // output in Fourier space
                    IGSTV.eta_pow_hat[p].data());
        }
        else {
            std::vector<double> eta_phy = ifft_real(eta_hat);
            for (size_t i = 0; i < eta_phy.size(); ++i) {
                IGSTV.eta_pow_hat[p][i] = std::complex<double>(std::pow(eta_phy[i] , p), 0.0);
            }
            IGSTV.eta_pow_hat[p] = fft(IGSTV.eta_pow_hat[p]);
        }
    }

    std::fill(IGSTV.beta_pow_hat[0].begin(), IGSTV.beta_pow_hat[0].end(), std::complex<double>(1.0, 0.0));
    if (DNO_Mb > 0) {
        IGSTV.beta_pow_hat[1] = beta_hat;
    }

    for (int p = 2; p <= DNO_Mb; ++p) {
        if (GIV.Dealiase_flag) {
            dealiaser(IGSTV.beta_pow_hat[p-1].data(), beta_hat.data(),
                    false,   // input in Fourier space
                    false,   // output in Fourier space
                    IGSTV.beta_pow_hat[p].data());
        }
        else {
            std::vector<double> beta_phy = ifft_real(beta_hat);
            for (size_t i = 0; i < beta_phy.size(); ++i) {
                IGSTV.beta_pow_hat[p][i] = std::complex<double>(std::pow(beta_phy[i] , p), 0.0);
            }
            IGSTV.beta_pow_hat[p] = fft(IGSTV.beta_pow_hat[p]);
        }
    }
    
}

void G1_MbM_2d(
            const std::vector<std::complex<double>>& eta_hat, 
            const std::vector<std::complex<double>>& beta_hat,
            const std::vector<std::complex<double>>& x_hat,
            std::vector<std::complex<double>>& G_result_hat,
            int DNO_M = GIV.DNO_M,
            int DNO_Mb = GIV.DNO_Mb);

void G2_MbM_2d(
            const std::vector<std::complex<double>>& eta_hat, 
            const std::vector<std::complex<double>>& beta_hat,
            const std::vector<std::complex<double>>& x_hat,
            std::vector<std::complex<double>>& G_result_hat,
            int DNO_M = GIV.DNO_M,
            int DNO_Mb = GIV.DNO_Mb);

void F1_F2_hat_2d(
            const std::vector<std::complex<double>>& eta_hat, 
            const std::vector<std::complex<double>>& beta_hat,
            const std::vector<std::complex<double>>& gamma_hat,
            const std::vector<std::complex<double>>& phis_hat,
            std::vector<std::complex<double>>& F1_hat,
            std::vector<std::complex<double>>& F2_hat);

#endif