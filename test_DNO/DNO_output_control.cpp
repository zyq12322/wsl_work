// DNO_output_control.cpp
#include "DNO_output_control.h"
#include "DNO_calc_OE.h"        // partial_eta_partial_t_2d, D_pow_hat, init_DNO_workspace_2d
#include "DNO_time_integral.h"  // DNOSystemParams 的完整定义
#include "var_op_function.h"    // fft_obj

#include <stdexcept>
#include <iomanip>

DNOOutputRecorder::DNOOutputRecorder(const DNOOutputConfig& cfg, double t0)
    : cfg_(cfg)
{
    record_ = cfg_.enabled
           && cfg_.interval > 0.0
           && !cfg_.x_axis.empty();

    if (!record_) {
        return;
    }

    if (!fft_obj) {
        throw std::runtime_error(
            "DNOOutputRecorder: fft_obj is not initialized");
    }

    eta_out_.open(cfg_.eta_file);
    xi_out_.open(cfg_.xi_file);
    beta_out_.open(cfg_.beta_file);
    EVI_out_.open(cfg_.EVI_file);
    if (!eta_out_ || !xi_out_ || !beta_out_ || !EVI_out_) {
        throw std::runtime_error(
            "DNOOutputRecorder: cannot open output files");
    }

    eta_out_  << std::setprecision(16);
    xi_out_   << std::setprecision(16);
    beta_out_ << std::setprecision(16);
    EVI_out_ << std::setprecision(16);

    // 表头: 0, x0, x1, ..., x_{N-1}
    eta_out_  << 0;
    xi_out_   << 0;
    beta_out_ << 0;
    for (double x : cfg_.x_axis) {
        eta_out_  << "," << x;
        xi_out_   << "," << x;
        beta_out_ << "," << x;
    }
    eta_out_  << "\n";
    xi_out_   << "\n";
    beta_out_ << "\n";

    next_record_t_ = t0;
}

void DNOOutputRecorder::record(
    double t,
    const std::vector<std::complex<double>>& eta_hat,
    const std::vector<std::complex<double>>& xi_hat,
    const std::vector<std::complex<double>>& beta_hat,
    const std::vector<std::complex<double>>& gamma_hat,
    const DNOSystemParams& p)
{
    if (!record_) return;

    const size_t N = cfg_.x_axis.size();
    if (eta_hat.size() != N || xi_hat.size() != N || beta_hat.size() != N || gamma_hat.size() != N) {
        throw std::invalid_argument(
            "DNOOutputRecorder::record: input size mismatch with x_axis");
    }
    if (!fft_obj || !dealiaser) {
        throw std::runtime_error(
            "DNOOutputRecorder::record: fft_obj or dealiaser is not initialized");
    }

    const std::complex<double> ImUnit(0.0, 1.0);

    // 频谱 -> 物理空间
    const auto eta_phy  = fft_obj->ifft(eta_hat);
    const auto xi_phy   = fft_obj->ifft(xi_hat);
    const auto beta_phy = fft_obj->ifft(beta_hat);

    std::vector<std::complex<double>> G_hat(N, {0.0, 0.0});

    init_DNO_workspace_2d(p.DNO_M, p.DNO_Mb, N, p.k_axis,
                          eta_hat, beta_hat, p.dealiase_flag);

    partial_eta_partial_t_2d(xi_hat, gamma_hat, G_hat,
                             p.h0, p.DNO_M, p.DNO_Mb, p.dealiase_flag);

    // --- 2) 逆变换 ---
    const std::vector<std::complex<double>> G_phy = fft_obj->ifft(G_hat);

    // --- 3) η_x 的谱微分 ---
    std::vector<std::complex<double>> eta_x_hat(N);
    for (size_t i = 0; i < N; ++i) {
        eta_x_hat[i] = ImUnit * D_pow_hat[1][i] * eta_hat[i];
    }
    const std::vector<std::complex<double>> eta_x_phy = fft_obj->ifft(eta_x_hat);

    // --- 4) 矩形法积分 ---
    //     非线性乘积 (ξ·G, η·η, ξ·η_x) 按 p.dealiase_flag 分情况处理：
    //       - dealiase_flag == true ：在谱空间通过 dealiaser 计算乘积，避免混叠
    //       - dealiase_flag == false：直接在物理空间逐点相乘
    const double dx = cfg_.x_axis[1] - cfg_.x_axis[0];
    double Ek = 0.0, Ep = 0.0, I = 0.0, V = 0.0;

    if (!p.dealiase_flag) {
        // ---- 物理空间直接逐点相乘 ----
        for (size_t i = 0; i < N; ++i) {
            const double eta_v  = eta_phy[i].real();
            const double beta_v  = beta_phy[i].real();
            const double xi_v   = xi_phy[i].real();
            const double G_v    = G_phy[i].real();
            const double etax_v = eta_x_phy[i].real();

            Ek += 0.5 * xi_v * G_v * dx;
            Ep += 0.5 * p.g * eta_v * eta_v * dx;
            I  +=       xi_v * etax_v * dx;
            V  += (eta_v - beta_v + p.h0) * dx;
        }
    } else {
        // ---- 谱空间通过 dealiaser 计算乘积，再逆变换积分 ----
        std::vector<std::complex<double>> xiG_hat(N);
        std::vector<std::complex<double>> eta2_hat(N);
        std::vector<std::complex<double>> xietax_hat(N);

        // ξ · (G ξ)
        (*dealiaser)(xi_hat.data(), G_hat.data(),
                     false, false, xiG_hat.data());
        // η · η
        (*dealiaser)(eta_hat.data(), eta_hat.data(),
                     false, false, eta2_hat.data());
        // ξ · η_x
        (*dealiaser)(xi_hat.data(), eta_x_hat.data(),
                     false, false, xietax_hat.data());

        const auto xiG_phy    = fft_obj->ifft(xiG_hat);
        const auto eta2_phy   = fft_obj->ifft(eta2_hat);
        const auto xietax_phy = fft_obj->ifft(xietax_hat);

        for (size_t i = 0; i < N; ++i) {
            Ek += 0.5 *       xiG_phy[i].real()    * dx;
            Ep += 0.5 * p.g * eta2_phy[i].real()   * dx;
            I  +=             xietax_phy[i].real() * dx;
            V  += (eta_phy[i].real() - beta_phy[i].real() + p.h0) * dx;
        }
    }

    eta_out_  << t;
    xi_out_   << t;
    beta_out_ << t;
    EVI_out_ << t;

    for (size_t i = 0; i < N; ++i) {
        eta_out_  << "," << eta_phy[i].real();
        xi_out_   << "," << xi_phy[i].real();
        beta_out_ << "," << beta_phy[i].real();
    }

    EVI_out_  << "," << Ek;
    EVI_out_  << "," << Ep;
    EVI_out_  << "," << I;
    EVI_out_  << "," << V;
    

    eta_out_  << "\n";
    xi_out_   << "\n";
    beta_out_ << "\n";
    EVI_out_ << "\n";

    eta_out_.flush();
    xi_out_.flush();
    beta_out_.flush();
    EVI_out_.flush();

    next_record_t_ = t + cfg_.interval;
}