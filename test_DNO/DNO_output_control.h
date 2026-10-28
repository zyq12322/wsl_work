// DNO_output_control.h
#ifndef DNO_OUTPUT_CONTROL_H
#define DNO_OUTPUT_CONTROL_H

#include <vector>
#include <complex>
#include <string>
#include <fstream>
#include <cstddef>

struct DNOSystemParams;

// ---------------------------------------------------------------------
// DNO 输出配置
//
// 每隔 interval 时间单位，把 η、ξ、β 的物理空间值写入文件。
// 文件格式：
//   第 1 行 : 0, x0, x1, ..., x_{N-1}
//   第 k 行 : t_k, f(x0), f(x1), ..., f(x_{N-1})
// ---------------------------------------------------------------------
struct DNOOutputConfig {
    bool                enabled   = false;
    double              interval  = 0.0;         // 记录时间间隔
    std::vector<double> x_axis;                  // 物理网格（长度 N）
    std::string         eta_file  = "eta.txt";
    std::string         xi_file   = "xi.txt";
    std::string         beta_file = "beta.txt";
    std::string         EVI_file = "EVI.txt";
};

// ---------------------------------------------------------------------
// DNO 输出记录器
//
// 使用方式：
//     DNOOutputRecorder rec(cfg, /*t0=*/0.0);          // 构造即建文件 + 写表头
//     rec.record(t0, eta_hat, xi_hat, beta_hat);       // 记录初始时刻（可选）
//     while (...) {
//         ... 积分一步 ...
//         if (rec.due(t)) {
//             std::vector<std::complex<double>> beta_hat(N), gamma_hat(N);
//             p.bc_func(t, beta_hat, gamma_hat);
//             rec.record(t, state.eta_hat, state.xi_hat, beta_hat);
//         }
//     }
//
// 说明：
//   - 若 cfg.enabled==false / interval<=0 / x_axis 为空，则记录器处于
//     disabled 状态，due() 恒返回 false，record() 不做任何事。
//   - 要求全局 fft_obj 已初始化（用于 ifft 到物理空间）。
//   - 三个输入向量必须与 x_axis 长度一致。
// ---------------------------------------------------------------------
class DNOOutputRecorder {
public:
    DNOOutputRecorder() = default;

    // 构造时若启用，则打开三个文件并写入表头行。
    // t0 用作首次记录时刻（next_record_t_ 初始化为 t0）。
    DNOOutputRecorder(const DNOOutputConfig& cfg, double t0);

    bool   enabled()   const noexcept { return record_; }

    // 是否应该在时刻 t 调用 record()（带 1e-12 容差）
    bool due(double t) const noexcept {
        return record_ && t >= next_record_t_ - 1e-12;
    }

    // 写出时刻 t 的一条记录，之后把 next_record_t_ 推进到 t + interval。
    // eta_hat / xi_hat / beta_hat / gamma_hat 均为傅里叶空间系数，长度 = N。
    void record(double t,
                const std::vector<std::complex<double>>& eta_hat,
                const std::vector<std::complex<double>>& xi_hat,
                const std::vector<std::complex<double>>& beta_hat,
                const std::vector<std::complex<double>>& gamma_hat,
                const DNOSystemParams& p);

private:
    DNOOutputConfig cfg_;
    bool            record_         = false;
    double          next_record_t_  = 0.0;

    std::ofstream   eta_out_;
    std::ofstream   xi_out_;
    std::ofstream   beta_out_;
    std::ofstream   EVI_out_;
};

#endif // DNO_OUTPUT_CONTROL_H