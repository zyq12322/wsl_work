#include "DNO_calc.h"
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <cassert>



void G1_m0_0_2d(// please confrim IGSTV.beta_pow_hat have been set correctly before call this function.
    const std::vector<std::complex<double>>& x_hat,//len == Nx(period), 2*Nx(wall)
    const std::vector<std::vector<std::complex<double>>>& G2_m0_m1_list,// shape:((DNO_Mb + 1)*(DNO_M+1), len(x_hat))
    std::vector<std::vector<std::complex<double>>>& G1_m0_m1_list,
    int DNO_M,
    int DNO_Mb)//fill this with G_1^{(m0,0)},m0=0,1,...,DNO_Mb
{
    size_t N = x_hat.size();

    // G1_m0_0_list[0, :] = x_hat * IGSV2D.k_axis * tanh(IGSV2D.k_axis * mean_depth)
    for (size_t i = 0; i < N; ++i) {
        G1_m0_m1_list[0][i] = x_hat[i] * IGSV2D.G00[i];
    }

    // m0_i = 1 ... DNO_Mb
    for (int m0_i = 1; m0_i <= DNO_Mb; ++m0_i) {
        // sum_term = zeros (长度 N)
        std::vector<std::complex<double>> sum_term(N, std::complex<double>(0.0, 0.0));

        // l 从 1 到 m0_i，步长为 2
        for (int l = 1; l <= m0_i; l += 2) {
            // odd_term = G2_m0_0_list[m0_i - l, :]
            const auto& row = G2_m0_m1_list[(m0_i - l) * (DNO_M + 1)];
            std::vector<std::complex<double>> odd_term(row.begin(), row.end());
            

            // odd_term = odd_term * IGSV2D.k_axis
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] *= IGSV2D.k_axis[i];
            }

            // odd_term /= factorial(l)
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] /= IGSV2D.fact_list[l];
            }

            // odd_term = a_product_bp(odd_term, beta_fft, l, False, False)
            if (GIV.Dealiase_flag) {
                (*dealiaser)(odd_term.data(), IGSTV.beta_pow_hat[l].data(),
                          false,   // input in Fourier space
                          false,   // output in Fourier space
                          odd_term.data());
            }
            else {
                odd_term = ifft(odd_term);
                std::vector<double> beta_phy = ifft_real(IGSTV.beta_pow_hat[1]);
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] *= std::pow(beta_phy[i] , l);
                }
                odd_term = fft(odd_term);
            }
            // odd_term *= (IGSV2D.k_axis ** l)
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] *= IGSV2D.D_pow_list[l][i];
            }

            // odd_term *= sech(mean_depth * IGSV2D.k_axis)
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] *= IGSV2D.sech_hD[i];
            }

            // sum_term -= odd_term
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] -= odd_term[i];
            }
        }

        // 将结果存入 G1_m0_0_list[m0_i, :]
        G1_m0_m1_list[m0_i * (DNO_M + 1)] = sum_term;
    }
}

void G2_m0_0_2d(// please confrim IGSTV.beta_pow_hat have been set correctly before call this function.
    const std::vector<std::complex<double>>& x_hat,//len == Nx(period), 2*Nx(wall)
    std::vector<std::vector<std::complex<double>>>& G2_m0_m1_list,
    int DNO_M,
    int DNO_Mb)//fill this with G_2^{(m0,0)},m0=0,1,...,DNO_Mb
{
    size_t N = x_hat.size();

    // G2_m0_0_list[0, :] = x_hat * sech(IGSV2D.k_axis * mean_depth)
    for (size_t i = 0; i < N; ++i) {
        G2_m0_m1_list[0][i] = x_hat[i] * IGSV2D.sech_hD[i];
    }

    // m0_i = 1 ... m0 
    for (int m0_i = 1; m0_i <= DNO_Mb; ++m0_i) {
        // sum_term = zeros (长度 N)
        std::vector<std::complex<double>> sum_term(N, std::complex<double>(0.0, 0.0));

        // l 从 1 到 m0_i，步长为 2
        for (int l = 1; l <= m0_i; l += 2) {
            // odd_term = G2_m0_0_list[m0_i - l, :]
            const auto& row = G2_m0_m1_list[(m0_i - l) * (DNO_M + 1)];
            std::vector<std::complex<double>> odd_term(row.begin(), row.end());

            // odd_term = odd_term * IGSV2D.k_axis
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] *= IGSV2D.k_axis[i];
            }

            // odd_term /= factorial(l)
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] /= IGSV2D.fact_list[l];
            }

            // odd_term = a_product_bp(odd_term, beta_fft, l, False, False)
            if (GIV.Dealiase_flag) {
                dealiaser(odd_term.data(), IGSTV.beta_pow_hat[l].data(),
                          false,   // input in Fourier space
                          false,   // output in Fourier space
                          odd_term.data());
            }
            else {
                odd_term = ifft(odd_term);
                std::vector<double> beta_phy = ifft_real(IGSTV.beta_pow_hat[1]);
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] *= std::pow(beta_phy[i] , l);
                }
                odd_term = fft(odd_term);
            }
            
            // odd_term *= (IGSV2D.k_axis ** (l-1))
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] *= IGSV2D.D_pow_list[l-1][i];
            }

            // odd_term *= tanh(mean_depth * IGSV2D.k_axis)
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] *= IGSV2D.tanh_hD[i];
            }

            // sum_term += odd_term
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] += odd_term[i];
            }
        }


        // l 从 2 到 m0_i，步长为 2
        for (int l = 2; l <= m0_i; l += 2) {
            // even_term = G2_m0_0_list[m0_i - l, :]
            const auto& row = G2_m0_m1_list[(m0_i - l) * (DNO_M + 1)];
            std::vector<std::complex<double>> even_term(row.begin(), row.end());

            // even_term = even_term * IGSV2D.k_axis
            for (size_t i = 0; i < N; ++i) {
                even_term[i] *= IGSV2D.k_axis[i];
            }

            // even_term /= factorial(l)
            for (size_t i = 0; i < N; ++i) {
                even_term[i] /= IGSV2D.fact_list[l];
            }

            // even_term = a_product_bp(even_term, beta_fft, l, False, False)
            if (GIV.Dealiase_flag) {
                dealiaser(even_term.data(), IGSTV.beta_pow_hat[l].data(),
                          false,   // input in Fourier space
                          false,   // output in Fourier space
                          even_term.data());
            }
            else {
                even_term = ifft(even_term);
                std::vector<double> beta_phy = ifft_real(IGSTV.beta_pow_hat[1]);
                for (size_t i = 0; i < N; ++i) {
                    even_term[i] *= std::pow(beta_phy[i] , l);
                }
                even_term = fft(even_term);
            }
            
            // even_term *= (IGSV2D.k_axis ** (l-1))
            for (size_t i = 0; i < N; ++i) {
                even_term[i] *= IGSV2D.D_pow_list[l-1][i];
            }

            // sum_term -= even_term
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] -= even_term[i];
            }
        }

        // 将结果存入 G2_m0_0_list[m0_i, :]
        G2_m0_m1_list[m0_i * (DNO_M + 1)] = sum_term;
    }
}

void G1_0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
    const std::vector<std::complex<double>>& x_hat,//len == Nx(period), 2*Nx(wall)
    std::vector<std::vector<std::complex<double>>>& G1_m0_m1_list,
    int DNO_M)//fill this with G_1^{(0,m1)},m1=0,1,...,DNO_M
{
    size_t N = x_hat.size();

    // G1_0_m1_list[0, :] = x_hat * IGSV2D.k_axis * tanh(IGSV2D.k_axis * mean_depth)
    for (size_t i = 0; i < N; ++i) {
        G1_m0_m1_list[0][i] = x_hat[i] * IGSV2D.G00[i];
    }

    // m1_i = 1 ... DNO_M 
    for (int m1_i = 1; m1_i <= DNO_M; ++m1_i) {
        // sum_term = zeros (长度 N)
        std::vector<std::complex<double>> sum_term(x_hat.begin(), x_hat.end());
        if (m1_i % 2 == 0) {
            // sum_term = sum_term * IGSV2D.k_axis
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] *= IGSV2D.k_axis[i];
            }
            // sum_term = a_product_bp(sum_term, eta_fft, m1_i, False, False)
            if (GIV.Dealiase_flag) {
                dealiaser(sum_term.data(), IGSTV.eta_pow_hat[m1_i].data(),
                        false,   // input in Fourier space
                        false,   // output in Fourier space
                        sum_term.data());
            }
            else {
                sum_term = ifft(sum_term);
                std::vector<double> eta_phy = ifft_real(IGSTV.eta_pow_hat[1]);
                for (size_t i = 0; i < N; ++i) {
                    sum_term[i] *= std::pow(eta_phy[i] , m1_i);
                }
                sum_term = fft(sum_term);
            }

            // sum_term = sum_term * IGSV2D.k_axis ** (m1_i-1) * G00
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] *= IGSV2D.D_pow_list[m1_i][i] * IGSV2D.tanh_hD[i];
            }
            // sum_term = sum_term / fact(m1_i)
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] /= IGSV2D.fact_list[m1_i];
            }

            // l 从 1 到 m1_i-1，步长为 2
            for (int l = 1; l <= m1_i-1; l += 2) {
                // odd_term = G1_0_m1_list[m1_i - l, :]
                const auto& row = G1_m0_m1_list[m1_i - l];
                std::vector<std::complex<double>> odd_term(row.begin(), row.end());

                // odd_term = a_product_bp(odd_term, eta_fft, l, False, False)
                if (GIV.Dealiase_flag) {
                    dealiaser(odd_term.data(), IGSTV.eta_pow_hat[l].data(),
                            false,   // input in Fourier space
                            false,   // output in Fourier space
                            odd_term.data());
                }
                else {
                    odd_term = ifft(odd_term);
                    std::vector<double> eta_phy = ifft_real(IGSTV.eta_pow_hat[1]);
                    for (size_t i = 0; i < N; ++i) {
                        odd_term[i] *= std::pow(eta_phy[i] , l);
                    }
                    odd_term = fft(odd_term);
                }

                // odd_term *= (IGSV2D.k_axis ** l) * tanh(h_0D)
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] *= IGSV2D.D_pow_list[l][i]  * IGSV2D.tanh_hD[i];
                }

                // odd_term /= factorial(l)
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] /= IGSV2D.fact_list[l];
                }

                // sum_term -= odd_term
                for (size_t i = 0; i < N; ++i) {
                    sum_term[i] -= odd_term[i];
                }
            }


            // l 从 2 到 m1_i，步长为 2
            for (int l = 2; l <= m1_i; l += 2) {
                // even_term = G1_m0_m1_list[m1_i - l, :]
                const auto& row = G1_m0_m1_list[m1_i - l];
                std::vector<std::complex<double>> even_term(row.begin(), row.end());

                // even_term = a_product_bp(even_term, eta_fft, l, False, False)
                if (GIV.Dealiase_flag) {
                    dealiaser(even_term.data(), IGSTV.eta_pow_hat[l].data(),
                            false,   // input in Fourier space
                            false,   // output in Fourier space
                            even_term.data());
                }
                else {
                    even_term = ifft(even_term);
                    std::vector<double> eta_phy = ifft_real(IGSTV.eta_pow_hat[1]);
                    for (size_t i = 0; i < N; ++i) {
                        even_term[i] *= std::pow(eta_phy[i] , l);
                    }
                    even_term = fft(even_term);
                }

                // even_term *= (IGSV2D.k_axis ** (l))
                for (size_t i = 0; i < N; ++i) {
                    even_term[i] *= IGSV2D.D_pow_list[l][i];
                }

                // even_term /= factorial(l)
                for (size_t i = 0; i < N; ++i) {
                    even_term[i] /= IGSV2D.fact_list[l];
                }

                // sum_term -= even_term
                for (size_t i = 0; i < N; ++i) {
                    sum_term[i] -= even_term[i];
                }
            }

            // 将结果存入 G1_0_m1_list[m1_i, :]
            G1_m0_m1_list[m1_i] = sum_term;
        }
        else {
            // sum_term = sum_term * IGSV2D.k_axis
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] *= IGSV2D.k_axis[i];
            }
            // sum_term = a_product_bp(sum_term, eta_fft, m1_i, False, False)
            if (GIV.Dealiase_flag) {
                dealiaser(sum_term.data(), IGSTV.eta_pow_hat[m1_i].data(),
                        false,   // input in Fourier space
                        false,   // output in Fourier space
                        sum_term.data());
            }
            else {
                sum_term = ifft(sum_term);
                std::vector<double> eta_phy = ifft_real(IGSTV.eta_pow_hat[1]);
                for (size_t i = 0; i < N; ++i) {
                    sum_term[i] *= std::pow(eta_phy[i] , m1_i);
                }
                sum_term = fft(sum_term);
            }

            // sum_term = sum_term * IGSV2D.k_axis ** (m1_i)
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] *= IGSV2D.D_pow_list[m1_i][i];
            }
            // sum_term = sum_term / fact(m1_i)
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] /= IGSV2D.fact_list[m1_i];
            }

            // l 从 1 到 m1_i-1，步长为 2
            for (int l = 1; l <= m1_i; l += 2) {
                // odd_term = G1_0_m1_list[m1_i - l, :]
                const auto& row = G1_m0_m1_list[m1_i - l];
                std::vector<std::complex<double>> odd_term(row.begin(), row.end());

                // odd_term = a_product_bp(odd_term, eta_fft, l, False, False)
                if (GIV.Dealiase_flag) {
                    dealiaser(odd_term.data(), IGSTV.eta_pow_hat[l].data(),
                            false,   // input in Fourier space
                            false,   // output in Fourier space
                            odd_term.data());
                }
                else {
                    odd_term = ifft(odd_term);
                    std::vector<double> eta_phy = ifft_real(IGSTV.eta_pow_hat[1]);
                    for (size_t i = 0; i < N; ++i) {
                        odd_term[i] *= std::pow(eta_phy[i] , l);
                    }
                    odd_term = fft(odd_term);
                }

                // odd_term *= (IGSV2D.k_axis ** l) * tanh(h_0D)
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] *= IGSV2D.D_pow_list[l][i]  * IGSV2D.tanh_hD[i];
                }

                // odd_term /= factorial(l)
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] /= IGSV2D.fact_list[l];
                }

                // sum_term -= odd_term
                for (size_t i = 0; i < N; ++i) {
                    sum_term[i] -= odd_term[i];
                }
            }


            // l 从 2 到 m1_i，步长为 2
            for (int l = 2; l <= m1_i-1; l += 2) {
                // even_term = G1_0_m1_list[m1_i - l, :]
                const auto& row = G1_m0_m1_list[m1_i - l];
                std::vector<std::complex<double>> even_term(row.begin(), row.end());

                // even_term = a_product_bp(even_term, eta_fft, l, False, False)
                if (GIV.Dealiase_flag) {
                    dealiaser(even_term.data(), IGSTV.eta_pow_hat[l].data(),
                            false,   // input in Fourier space
                            false,   // output in Fourier space
                            even_term.data());
                }
                else {
                    even_term = ifft(even_term);
                    std::vector<double> eta_phy = ifft_real(IGSTV.eta_pow_hat[1]);
                    for (size_t i = 0; i < N; ++i) {
                        even_term[i] *= std::pow(eta_phy[i] , l);
                    }
                    even_term = fft(even_term);
                }

                // even_term *= (IGSV2D.k_axis ** (l))
                for (size_t i = 0; i < N; ++i) {
                    even_term[i] *= IGSV2D.D_pow_list[l][i];
                }

                // even_term /= factorial(l)
                for (size_t i = 0; i < N; ++i) {
                    even_term[i] /= IGSV2D.fact_list[l];
                }

                // sum_term -= even_term
                for (size_t i = 0; i < N; ++i) {
                    sum_term[i] -= even_term[i];
                }
            }

            // 将结果存入 G1_m0_m1_list[m1_i, :]
            G1_m0_m1_list[m1_i] = sum_term;
        }
    }
}

void G2_0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
    const std::vector<std::complex<double>>& x_hat,//len == Nx(period), 2*Nx(wall)
    std::vector<std::vector<std::complex<double>>>& G1_m0_m1_list,
    std::vector<std::vector<std::complex<double>>>& G2_m0_m1_list,
    int DNO_M)//fill this with G_2^{(0,m1)},m1=0,1,...,DNO_M
{
    size_t N = x_hat.size();

    // G2_0_m1_list[0, :] = x_hat * sech(IGSV2D.k_axis * mean_depth)
    for (size_t i = 0; i < N; ++i) {
        G2_m0_m1_list[0][i] = x_hat[i] * IGSV2D.sech_hD[i];
    }

    // m1_i = 1 ... DNO_M
    for (int m1_i = 1; m1_i <= DNO_M; ++m1_i) {
        if (m1_i % 2 == 0){
            // sum_term = x_hat (长度 N)
            std::vector<std::complex<double>> sum_term(x_hat.begin(), x_hat.end());
            // sum_term = sum_term * IGSV2D.k_axis
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] *= IGSV2D.k_axis[i];
            }
            // sum_term = a_product_bp(sum_term, eta_fft, m1_i, False, False)
            if (GIV.Dealiase_flag) {
                dealiaser(sum_term.data(), IGSTV.eta_pow_hat[m1_i].data(),
                        false,   // input in Fourier space
                        false,   // output in Fourier space
                        sum_term.data());
            }
            else {
                sum_term = ifft(sum_term);
                std::vector<double> eta_phy = ifft_real(IGSTV.eta_pow_hat[1]);
                for (size_t i = 0; i < N; ++i) {
                    sum_term[i] *= std::pow(eta_phy[i] , m1_i);
                }
                sum_term = fft(sum_term);
            }

            // sum_term = sum_term * IGSV2D.k_axis ** (m1_i-1) * sech(h_0D)
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] *= IGSV2D.D_pow_list[m1_i - 1][i] * IGSV2D.sech_hD[i];
            }

            for (size_t i = 0; i < N; ++i) {
                sum_term[i] /= IGSV2D.fact_list[m1_i];
            }

            // l 从 1 到 m0_i，步长为 2
            for (int l = 1; l <= m1_i-1; l += 2) {
                // odd_term = G1_0_m1_list[m1_i - l, :]
                const auto& row = G1_m0_m1_list[m1_i - l];
                std::vector<std::complex<double>> odd_term(row.begin(), row.end());

                // odd_term = a_product_bp(odd_term, eta_fft, l, False, False)
                if (GIV.Dealiase_flag) {
                    dealiaser(odd_term.data(), IGSTV.eta_pow_hat[l].data(),
                            false,   // input in Fourier space
                            false,   // output in Fourier space
                            odd_term.data());
                }
                else {
                    odd_term = ifft(odd_term);
                    std::vector<double> eta_phy = ifft_real(IGSTV.eta_pow_hat[1]);
                    for (size_t i = 0; i < N; ++i) {
                        odd_term[i] *= std::pow(eta_phy[i] , l);
                    }
                    odd_term = fft(odd_term);
                }

                // odd_term *= (IGSV2D.k_axis ** (l-1)) * sech(h_0D)
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] *= IGSV2D.D_pow_list[l-1][i] * IGSV2D.sech_hD[i];
                }
                // odd_term /= factorial(l)
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] /= IGSV2D.fact_list[l];
                }

                // sum_term -= odd_term
                for (size_t i = 0; i < N; ++i) {
                    sum_term[i] -= odd_term[i];
                }
            }
            G2_m0_m1_list[m1_i] = sum_term;
        }
        else {
            // sum_term = zeros (长度 N)
            std::vector<std::complex<double>> sum_term(N, std::complex<double>(0.0, 0.0));
            
            // l 从 1 到 m0_i，步长为 2
            for (int l = 1; l <= m1_i; l += 2) {
                // odd_term = G1_0_m1_list[m1_i - l, :]
                const auto& row = G1_m0_m1_list[m1_i - l];
                std::vector<std::complex<double>> odd_term(row.begin(), row.end());

                // odd_term = a_product_bp(odd_term, eta_fft, l, False, False)
                if (GIV.Dealiase_flag) {
                    dealiaser(odd_term.data(), IGSTV.eta_pow_hat[l].data(),
                            false,   // input in Fourier space
                            false,   // output in Fourier space
                            odd_term.data());
                }
                else {
                    odd_term = ifft(odd_term);
                    std::vector<double> eta_phy = ifft_real(IGSTV.eta_pow_hat[1]);
                    for (size_t i = 0; i < N; ++i) {
                        odd_term[i] *= std::pow(eta_phy[i] , l);
                    }
                    odd_term = fft(odd_term);
                }

                // odd_term *= (IGSV2D.k_axis ** (l-1)) * sech(h_0D)
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] *= IGSV2D.D_pow_list[l-1][i] * IGSV2D.sech_hD[i];
                }
                // odd_term /= factorial(l)
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] /= IGSV2D.fact_list[l];
                }

                // sum_term -= odd_term
                for (size_t i = 0; i < N; ++i) {
                    sum_term[i] -= odd_term[i];
                }
            }
            G2_m0_m1_list[m1_i] = sum_term;
        }
    
    }
}

void G1_m0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
    int m0,
    const std::vector<std::complex<double>>& x_hat,//len == Nx(period), 2*Nx(wall)
    std::vector<std::vector<std::complex<double>>>& G2_m0_m1_list,
    std::vector<std::vector<std::complex<double>>>& G1_m0_m1_list,
    int DNO_M,
    int DNO_Mb)//fill this with G_1^{(m0,m1)},m1=1,...,DNO_M
{
    if (m0 <= 0) {
        throw std::invalid_argument("value must be positive in G1_m0_m1_2d");
    }
    size_t N = x_hat.size();

    // m1_i = 1 ... DNO_M 
    for (int m1_i = 1; m1_i <= DNO_M; ++m1_i) {
        // sum_term = zeros (长度 N)
        std::vector<std::complex<double>> sum_term(N, std::complex<double>(0.0, 0.0));

        // l 从 1 到 m1_i，步长为 2
        for (int l = 1; l <= m1_i; l += 2) {
            // odd_term = G1_m0_m1_list[m0, m1_i - l, :]
            const auto& row = G1_m0_m1_list[m0*(DNO_M + 1) + m1_i - l];
            std::vector<std::complex<double>> odd_term(row.begin(), row.end());

            // odd_term = a_product_bp(odd_term, eta_fft, l, False, False)
            if (GIV.Dealiase_flag) {
                dealiaser(odd_term.data(), IGSTV.eta_pow_hat[l].data(),
                        false,   // input in Fourier space
                        false,   // output in Fourier space
                        odd_term.data());
            }
            else {
                odd_term = ifft(odd_term);
                std::vector<double> eta_phy = ifft_real(IGSTV.eta_pow_hat[1]);
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] *= std::pow(eta_phy[i] , l);
                }
                odd_term = fft(odd_term);
            }

            // odd_term *= (IGSV2D.k_axis ** l) * tanh(h_0D)
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] *= IGSV2D.D_pow_list[l][i]  * IGSV2D.tanh_hD[i];
            }

            // odd_term /= factorial(l)
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] /= IGSV2D.fact_list[l];
            }

            // sum_term -= odd_term
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] -= odd_term[i];
            }
        }


        // l 从 2 到 m1_i，步长为 2
        for (int l = 2; l <= m1_i; l += 2) {
            // even_term = G1_m0_m1_list[m0, m1_i - l, :]
            const auto& row = G1_m0_m1_list[m0 *(DNO_M + 1) + (m1_i - l)];
            std::vector<std::complex<double>> even_term(row.begin(), row.end());

            // even_term = a_product_bp(even_term, eta_fft, l, False, False)
            if (GIV.Dealiase_flag) {
                dealiaser(even_term.data(), IGSTV.eta_pow_hat[l].data(),
                        false,   // input in Fourier space
                        false,   // output in Fourier space
                        even_term.data());
            }
            else {
                even_term = ifft(even_term);
                std::vector<double> eta_phy = ifft_real(IGSTV.eta_pow_hat[1]);
                for (size_t i = 0; i < N; ++i) {
                    even_term[i] *= std::pow(eta_phy[i] , l);
                }
                even_term = fft(even_term);
            }

            // even_term *= (IGSV2D.k_axis ** (l))
            for (size_t i = 0; i < N; ++i) {
                even_term[i] *= IGSV2D.D_pow_list[l][i];
            }

            // even_term /= factorial(l)
            for (size_t i = 0; i < N; ++i) {
                even_term[i] /= IGSV2D.fact_list[l];
            }

            // sum_term -= even_term
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] -= even_term[i];
            }
        }

        // l 从 1 到 m0，步长为 2
        for (int l = 1; l <= m0; l += 2) {
            // odd_term = G2_m0_m1_list[m0-l, m1_i, :]
            const auto& row = G2_m0_m1_list[(m0-l) * (DNO_M + 1) + m1_i];
            std::vector<std::complex<double>> odd_term(row.begin(), row.end());

            // odd_term *= IGSV2D.k_axis
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] *= IGSV2D.k_axis[i];
            }

            // odd_term /= factorial(l)
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] /= IGSV2D.fact_list[l];
            }

            // odd_term = a_product_bp(odd_term, beta_fft, l, False, False)
            if (GIV.Dealiase_flag) {
                dealiaser(odd_term.data(), IGSTV.beta_pow_hat[l].data(),
                        false,   // input in Fourier space
                        false,   // output in Fourier space
                        odd_term.data());
            }
            else {
                odd_term = ifft(odd_term);
                std::vector<double> beta_phy = ifft_real(IGSTV.beta_pow_hat[1]);
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] *= std::pow(beta_phy[i] , l);
                }
                odd_term = fft(odd_term);
            }

            // odd_term *= (IGSV2D.k_axis ** l) * sech(h_0D)
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] *= IGSV2D.D_pow_list[l][i]  * IGSV2D.sech_hD[i];
            }

            // sum_term -= odd_term
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] -= odd_term[i];
            }
        }

        // 将结果存入 G1_m0_m1_list[m0, m1_i, :]
        G1_m0_m1_list[m0*(DNO_M + 1) + m1_i] = sum_term;
    }
}

void G2_m0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
    int m0,
    const std::vector<std::complex<double>>& x_hat,//len == Nx(period), 2*Nx(wall)
    std::vector<std::vector<std::complex<double>>>& G1_m0_m1_list,
    std::vector<std::vector<std::complex<double>>>& G2_m0_m1_list,
    int DNO_M,
    int DNO_Mb)//fill this with G_2^{(m0,m1)},m1=1,...,DNO_M
{
    if (m0 <= 0) {
        throw std::invalid_argument("value must be positive in G2_m0_m1_2d");
    }
    size_t N = x_hat.size();

    // m1_i = 1 ... DNO_M 
    for (int m1_i = 1; m1_i <= DNO_M; ++m1_i) {
        // sum_term = zeros (长度 N)
        std::vector<std::complex<double>> sum_term(N, std::complex<double>(0.0, 0.0));

        // l 从 1 到 m1_i-1，步长为 2
        for (int l = 1; l <= m1_i; l += 2) {
            // odd_term = G1_m0_m1_list[m0, m1_i - l, :]
            const auto& row = G1_m0_m1_list[m0*(DNO_M + 1) + m1_i - l];
            std::vector<std::complex<double>> odd_term(row.begin(), row.end());

            // odd_term = a_product_bp(odd_term, eta_fft, l, False, False)
            if (GIV.Dealiase_flag) {
                dealiaser(odd_term.data(), IGSTV.eta_pow_hat[l].data(),
                        false,   // input in Fourier space
                        false,   // output in Fourier space
                        odd_term.data());
            }
            else {
                odd_term = ifft(odd_term);
                std::vector<double> eta_phy = ifft_real(IGSTV.eta_pow_hat[1]);
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] *= std::pow(eta_phy[i] , l);
                }
                odd_term = fft(odd_term);
            }

            // odd_term *= (IGSV2D.k_axis ** l) * tanh(h_0D)
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] *= IGSV2D.D_pow_list[l-1][i]  * IGSV2D.sech_hD[i];
            }

            // odd_term /= factorial(l)
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] /= IGSV2D.fact_list[l];
            }

            // sum_term -= odd_term
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] -= odd_term[i];
            }
        }


        // l 从 2 到 m1_i，步长为 2
        for (int l = 2; l <= m0; l += 2) {
            // even_term = G2_m0_m1_list[m0-l, m1_i, :]
            const auto& row = G2_m0_m1_list[(m0-l) *(DNO_M + 1) + (m1_i)];
            std::vector<std::complex<double>> even_term(row.begin(), row.end());

            // even_term *= (IGSV2D.k_axis)
            for (size_t i = 0; i < N; ++i) {
                even_term[i] *= IGSV2D.k_axis[i];
            }

            // even_term /= factorial(l)
            for (size_t i = 0; i < N; ++i) {
                even_term[i] /= IGSV2D.fact_list[l];
            }

            // even_term = a_product_bp(even_term, beta_fft, l, False, False)
            if (GIV.Dealiase_flag) {
                dealiaser(even_term.data(), IGSTV.beta_pow_hat[l].data(),
                        false,   // input in Fourier space
                        false,   // output in Fourier space
                        even_term.data());
            }
            else {
                even_term = ifft(even_term);
                std::vector<double> beta_phy = ifft_real(IGSTV.beta_pow_hat[1]);
                for (size_t i = 0; i < N; ++i) {
                    even_term[i] *= std::pow(beta_phy[i] , l);
                }
                even_term = fft(even_term);
            }

            // even_term *= (IGSV2D.k_axis ** (l-1))
            for (size_t i = 0; i < N; ++i) {
                even_term[i] *= IGSV2D.D_pow_list[l-1][i];
            }

            // sum_term -= even_term
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] -= even_term[i];
            }
        }

        // l 从 1 到 m0，步长为 2
        for (int l = 1; l <= m0; l += 2) {
            // odd_term = G2_m0_m1_list[m0-l, m1_i, :]
            const auto& row = G2_m0_m1_list[(m0-l) * (DNO_M + 1) + m1_i];
            std::vector<std::complex<double>> odd_term(row.begin(), row.end());

            // odd_term *= IGSV2D.k_axis
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] *= IGSV2D.k_axis[i];
            }

            // odd_term /= factorial(l)
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] /= IGSV2D.fact_list[l];
            }

            // odd_term = a_product_bp(odd_term, beta_fft, l, False, False)
            if (GIV.Dealiase_flag) {
                dealiaser(odd_term.data(), IGSTV.beta_pow_hat[l].data(),
                        false,   // input in Fourier space
                        false,   // output in Fourier space
                        odd_term.data());
            }
            else {
                odd_term = ifft(odd_term);
                std::vector<double> beta_phy = ifft_real(IGSTV.beta_pow_hat[1]);
                for (size_t i = 0; i < N; ++i) {
                    odd_term[i] *= std::pow(beta_phy[i] , l);
                }
                odd_term = fft(odd_term);
            }

            // odd_term *= (IGSV2D.k_axis ** (l-1)) * tanh(h_0D)
            for (size_t i = 0; i < N; ++i) {
                odd_term[i] *= IGSV2D.D_pow_list[l-1][i]  * IGSV2D.tanh_hD[i];
            }

            // sum_term -= odd_term
            for (size_t i = 0; i < N; ++i) {
                sum_term[i] += odd_term[i];
            }
        }

        // 将结果存入 G2_m0_m1_list[m0, m1_i, :]
        G2_m0_m1_list[m0*(DNO_M + 1) + m1_i] = sum_term;
    }
}

void G1_MbM_2d(
            const std::vector<std::complex<double>>& eta_hat, 
            const std::vector<std::complex<double>>& beta_hat,
            const std::vector<std::complex<double>>& x_hat,
            std::vector<std::complex<double>>& G_result_hat,
            int DNO_M,
            int DNO_Mb) 
{
    size_t N = x_hat.size();
    if (eta_hat.size() != N || beta_hat.size() != N || G_result_hat.size() != N) {
        throw std::runtime_error("Size mismatch: eta_hat, beta_hat, and G_result_hat must all have size N");
    }
    std::vector<std::vector<std::complex<double>>> G1_m0_m1_list;
    std::vector<std::vector<std::complex<double>>> G2_m0_m1_list;
    resize_G_m0_m1_list_2d(G1_m0_m1_list, DNO_Mb, DNO_M, N);
    resize_G_m0_m1_list_2d(G2_m0_m1_list, DNO_Mb, DNO_M, N);
    resize_eta_beta_pow_hat(DNO_M, DNO_Mb, N);

    std::vector<double> beta_phy = ifft_real(beta_hat);
    double max_abs_beta = 0.0;
    for (double val : beta_phy) {
        double av = std::fabs(val);
        if (av > max_abs_beta) max_abs_beta = av;
    }
    int beta_scale = static_cast<int>(std::round(max_abs_beta));
    if (beta_scale == 0) beta_scale = 1;
    const std::complex<double> scale_factor(1.0 / static_cast<double>(beta_scale), 0.0);
    std::vector<std::complex<double>> beta_hat_scaled = beta_hat;
    for (size_t i = 0; i < beta_hat_scaled.size(); ++i) {
        beta_hat_scaled[i] *= scale_factor;
    }

    set_eta_beta_pow_hat(dealiaser, 
                         eta_hat, 
                         beta_hat_scaled,
                         DNO_Mb,
                         DNO_M);
    G2_m0_0_2d(// please confrim IGSTV.beta_pow_hat have been set correctly before call this function.
                dealiaser,
                x_hat,//len == Nx(period), 2*Nx(wall)
                G2_m0_m1_list,
                DNO_M,
                DNO_Mb);
    G1_m0_0_2d(
                dealiaser,
                x_hat,//len == Nx(period), 2*Nx(wall)
                G2_m0_m1_list,
                G1_m0_m1_list,
                DNO_M,
                DNO_Mb);
    G1_0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
                dealiaser,
                x_hat,//len == Nx(period), 2*Nx(wall)
                G1_m0_m1_list,
                DNO_M);//fill this with G_1^{(0,m1)},m1=0,1,...,DNO_M
    G2_0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
                dealiaser,
                x_hat,//len == Nx(period), 2*Nx(wall)
                G1_m0_m1_list,
                G2_m0_m1_list,
                DNO_M);

    for (int m0_i = 1;m0_i <= DNO_Mb; ++m0_i ) {
        G1_m0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
                    dealiaser,
                    m0_i,
                    x_hat,//len == Nx(period), 2*Nx(wall)
                    G2_m0_m1_list,
                    G1_m0_m1_list,
                    DNO_M,
                    DNO_Mb);//fill this with G_1^{(m0,m1)},m1=1,...,DNO_M
        G2_m0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
                    dealiaser,
                    m0_i,
                    x_hat,//len == Nx(period), 2*Nx(wall)
                    G1_m0_m1_list,
                    G2_m0_m1_list,
                    DNO_M,
                    DNO_Mb);
    }
    std::fill(G_result_hat.begin(), G_result_hat.end(), std::complex<double>(0.0, 0.0));
    
    for (int mi = 0; mi < (DNO_M + 1)* (DNO_Mb + 1); ++mi) {
        int m0_i = mi / (DNO_M + 1);
        for (size_t ind = 0 ; ind < G_result_hat.size(); ++ind) {
            G_result_hat[ind] += G1_m0_m1_list[mi][ind] * std::pow(static_cast<double>(beta_scale), m0_i);
        }
    }

}

void G2_MbM_2d(
            const std::vector<std::complex<double>>& eta_hat, 
            const std::vector<std::complex<double>>& beta_hat,
            const std::vector<std::complex<double>>& x_hat,
            std::vector<std::complex<double>>& G_result_hat,
            int DNO_M,
            int DNO_Mb) 
{
    size_t N = x_hat.size();
    if (eta_hat.size() != N || beta_hat.size() != N || G_result_hat.size() != N) {
        throw std::runtime_error("Size mismatch: eta_hat, beta_hat, and G_result_hat must all have size N");
    }
    std::vector<std::vector<std::complex<double>>> G1_m0_m1_list;
    std::vector<std::vector<std::complex<double>>> G2_m0_m1_list;
    resize_G_m0_m1_list_2d(G1_m0_m1_list, DNO_Mb, DNO_M, N);
    resize_G_m0_m1_list_2d(G2_m0_m1_list, DNO_Mb, DNO_M, N);
    resize_eta_beta_pow_hat(DNO_M, DNO_Mb, N);

    std::vector<double> beta_phy = ifft_real(beta_hat);
    double max_abs_beta = 0.0;
    for (double val : beta_phy) {
        double av = std::fabs(val);
        if (av > max_abs_beta) max_abs_beta = av;
    }
    int beta_scale = static_cast<int>(std::round(max_abs_beta));
    if (beta_scale == 0) beta_scale = 1;
    const std::complex<double> scale_factor(1.0 / static_cast<double>(beta_scale), 0.0);
    std::vector<std::complex<double>> beta_hat_scaled = beta_hat;
    for (size_t i = 0; i < beta_hat_scaled.size(); ++i) {
        beta_hat_scaled[i] *= scale_factor;
    }

    set_eta_beta_pow_hat(dealiaser, 
                         eta_hat, 
                         beta_hat_scaled,
                         DNO_Mb,
                         DNO_M);
    G2_m0_0_2d(// please confrim IGSTV.beta_pow_hat have been set correctly before call this function.
                dealiaser,
                x_hat,//len == Nx(period), 2*Nx(wall)
                G2_m0_m1_list,
                DNO_M,
                DNO_Mb);
    G1_m0_0_2d(
                dealiaser,
                x_hat,//len == Nx(period), 2*Nx(wall)
                G2_m0_m1_list,
                G1_m0_m1_list,
                DNO_M,
                DNO_Mb);
    G1_0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
                dealiaser,
                x_hat,//len == Nx(period), 2*Nx(wall)
                G1_m0_m1_list,
                DNO_M);//fill this with G_1^{(0,m1)},m1=0,1,...,DNO_M
    G2_0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
                dealiaser,
                x_hat,//len == Nx(period), 2*Nx(wall)
                G1_m0_m1_list,
                G2_m0_m1_list,
                DNO_M);

    for (int m0_i = 1;m0_i <= DNO_Mb; ++m0_i ) {
        G1_m0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
                    dealiaser,
                    m0_i,
                    x_hat,//len == Nx(period), 2*Nx(wall)
                    G2_m0_m1_list,
                    G1_m0_m1_list,
                    DNO_M,
                    DNO_Mb);//fill this with G_1^{(m0,m1)},m1=1,...,DNO_M
        G2_m0_m1_2d(// please confrim IGSTV.eta_pow_hat have been set correctly before call this function.
                    dealiaser,
                    m0_i,
                    x_hat,//len == Nx(period), 2*Nx(wall)
                    G1_m0_m1_list,
                    G2_m0_m1_list,
                    DNO_M,
                    DNO_Mb);
    }
    std::fill(G_result_hat.begin(), G_result_hat.end(), std::complex<double>(0.0, 0.0));
    for (int mi = 0; mi < (DNO_M + 1)* (DNO_Mb + 1); ++mi) {
        int m0_i = mi / (DNO_M + 1);
        for (size_t ind = 0 ; ind < G_result_hat.size(); ++ind) {
            G_result_hat[ind] += G2_m0_m1_list[mi][ind] * std::pow(static_cast<double>(beta_scale), m0_i);
        }
    }
}

void F1_F2_hat_2d(
            const std::vector<std::complex<double>>& eta_hat, 
            const std::vector<std::complex<double>>& beta_hat,
            const std::vector<std::complex<double>>& gamma_hat,
            const std::vector<std::complex<double>>& phis_hat,
            std::vector<std::complex<double>>& F1_hat,
            std::vector<std::complex<double>>& F2_hat)
{
    size_t N = phis_hat.size();
    if (eta_hat.size() != N || beta_hat.size() != N || gamma_hat.size() != N) {
        throw std::runtime_error("Size mismatch: eta_hat, beta_hat, and gamma_hat must all have size N");
    }
    if (F1_hat.size() != N || F2_hat.size() != N) {
        throw std::runtime_error("Size mismatch: phis_hat, F1_hat and F2_hat must all have size N");
    }

    std::vector<std::complex<double>> G2_gamma_hat(N);
    std::vector<std::complex<double>> G_hat(N);
    G1_MbM_2d(dealiaser, 
                eta_hat, 
                beta_hat,
                phis_hat,
                G_hat,
                GIV.DNO_M,
                GIV.DNO_Mb) ;
    if (VB.time_dependency){
        G2_MbM_2d(dealiaser, 
                    eta_hat, 
                    beta_hat,
                    gamma_hat,
                    G2_gamma_hat) ;
        for (size_t i = 0; i< N; ++i) {
            G_hat[i] += G2_gamma_hat[i];
        }
    }
    double adj = 1.0;
    if (GIV.NL_Ta > 1e-3) {
        adj = 1.0 - std::exp( - std::pow( std::abs(IGSTV.t_count) / GIV.NL_Ta, 4 ) );
    }
    std::vector<std::complex<double>> eta_x_hat(N);
    std::vector<std::complex<double>> phis_x_hat(N);
    std::complex<double> imag(0.0, 1.0);
    for (size_t ind = 0 ; ind < N; ++ind) {
        F1_hat[ind] = (G_hat[ind] - phis_hat[ind] * IGSV2D.G00[ind]) * adj;
        eta_x_hat[ind] = eta_hat[ind] * imag * IGSV2D.k_axis[ind];
        phis_x_hat[ind] = phis_hat[ind] * imag * IGSV2D.k_axis[ind];
    }

    std::vector<double> eta_x_phy(N);
    std::vector<double> phis_x_phy(N);
    std::vector<std::complex<double>> G_phy(N);
    eta_x_phy = ifft_real(eta_x_hat);
    phis_x_phy = ifft_real(phis_x_hat);
    G_phy = ifft(G_hat);

    for (size_t ind = 0 ; ind < N; ++ind) {
        F2_hat[ind] = (-0.5 / ( 1+ eta_x_phy[ind] * eta_x_phy[ind] ) 
        * ( phis_x_phy[ind] * phis_x_phy[ind] - G_phy[ind] * G_phy[ind]
        - 2.0 * eta_x_phy[ind] * phis_x_phy[ind] * G_phy[ind] ) ) * adj;
    }

    F2_hat = fft(F2_hat);
}
