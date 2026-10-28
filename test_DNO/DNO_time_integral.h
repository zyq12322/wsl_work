#ifndef DNO_TIME_INTEGRAL_H
#define DNO_TIME_INTEGRAL_H

#include "DNO_calc_OE.h"
#include "DNO_output_control.h"

#include <vector>
#include <complex>
#include <functional>
#include <cstddef>

// ---------------------------------------------------------------------
// β、γ 时间回调：给定时刻 t，填充 beta_hat 和 gamma_hat
//
// 说明：
//   - 若 β、γ 为常数，回调里直接填入预置值即可；
//   - 若 β 随时间变化且 γ = ∂β/∂t 有解析表达式，可在同一回调中计算；
//   - 若 γ 由外部数值方法提供，也在这里统一返回。
// ---------------------------------------------------------------------
using DNOBcFunc = std::function<void(
    double t,
    std::vector<std::complex<double>>& beta_hat,
    std::vector<std::complex<double>>& gamma_hat)>;

// ---------------------------------------------------------------------
// DNO 状态：η、ξ 的傅里叶系数
// ---------------------------------------------------------------------
struct DNOState {
    std::vector<std::complex<double>> eta_hat;
    std::vector<std::complex<double>> xi_hat;
};

// ---------------------------------------------------------------------
// 系统参数：物理常数、阶数、波数轴、β/γ 回调
//
// 调用者必须提前：
//   1. 创建全局 fft_obj / dealiaser；
//   2. 调用 resize_DNO_workspace(DNO_M, DNO_Mb, N) 一次。
// 本模块不会主动 resize。
// ---------------------------------------------------------------------
struct DNOSystemParams {
    double h0;                      // 平均水深
    double g;                       // 重力加速度
    int    DNO_M;                   // η 最大幂次
    int    DNO_Mb;                  // β 最大幂次
    bool   dealiase_flag;           // 是否使用反混叠
    std::vector<double> k_axis;     // 角波数轴 (rad/m)，长度 N
    DNOBcFunc bc_func;              // β、γ 时间回调
};

// ---------------------------------------------------------------------
// 单步尝试结果
// ---------------------------------------------------------------------
struct RK45AttemptResult {
    bool   accepted;    // 该步是否被接受
    double new_dt;      // 建议的下一步步长
    double err_norm;    // 加权 RMS 误差范数
};

// ---------------------------------------------------------------------
// Dormand-Prince RK45(5(4)) 单步
//
//   - 接受：state 被更新为 t+dt 处的 5 阶解；
//   - 拒绝：state 保持不变，调用者需以更小的 dt 重试。
//
// 参数：
//   t, dt            当前时刻与尝试步长
//   p                系统参数（含 β/γ 回调）
//   state            输入/输出：当前状态
//   rtol, atol       相对/绝对误差容限
//   safety           步长安全因子（默认 0.9）
//   min_factor       单步最小步长缩放（默认 0.2）
//   max_factor       单步最大步长缩放（默认 5.0）
// ---------------------------------------------------------------------
RK45AttemptResult rk45_step(
    double t, double dt,
    const DNOSystemParams& p,
    DNOState& state,
    double rtol, double atol,
    double safety     = 0.9,
    double min_factor = 0.2,
    double max_factor = 5.0);

// ---------------------------------------------------------------------
// RK45 阶段缓冲池：只在首次调用时分配，之后原地复用。
// ---------------------------------------------------------------------
struct RK45BufferCache
{
    std::vector<std::complex<double>> k1e, k1x;
    std::vector<std::complex<double>> k2e, k2x;
    std::vector<std::complex<double>> k3e, k3x;
    std::vector<std::complex<double>> k4e, k4x;
    std::vector<std::complex<double>> k5e, k5x;
    std::vector<std::complex<double>> k6e, k6x;
    std::vector<std::complex<double>> k7e, k7x;
    std::vector<std::complex<double>> eta_tmp, xi_tmp;
    std::vector<std::complex<double>> eta_new, xi_new;
    std::vector<std::complex<double>> eta_err, xi_err;

    // ---- FSAL 缓存状态 ----
    double cached_k1_t  = 0.0;      // 缓存 k1 所在的时刻
    bool   has_cached_k1 = false;   // 缓存是否有效

    DNOBufferCache dno;

    void ensure_size(size_t N, int DNO_M, int DNO_Mb)
    {
        auto needs = [N](const std::vector<std::complex<double>>& v) {
            return v.size() != N;
        };
        if (needs(k1e) || needs(k1x) ||
            needs(k2e) || needs(k2x) ||
            needs(k3e) || needs(k3x) ||
            needs(k4e) || needs(k4x) ||
            needs(k5e) || needs(k5x) ||
            needs(k6e) || needs(k6x) ||
            needs(k7e) || needs(k7x) ||
            needs(eta_tmp) || needs(xi_tmp) ||
            needs(eta_new) || needs(xi_new) ||
            needs(eta_err) || needs(xi_err))
        {
            auto alloc_N = [N](std::vector<std::complex<double>>& v){
                v.resize(N);
            };
            alloc_N(k1e); alloc_N(k1x); alloc_N(k2e); alloc_N(k2x);
            alloc_N(k3e); alloc_N(k3x); alloc_N(k4e); alloc_N(k4x);
            alloc_N(k5e); alloc_N(k5x); alloc_N(k6e); alloc_N(k6x);
            alloc_N(k7e); alloc_N(k7x);
            alloc_N(eta_tmp); alloc_N(xi_tmp);
            alloc_N(eta_new); alloc_N(xi_new);
            alloc_N(eta_err); alloc_N(xi_err);

            has_cached_k1 = false;   // ← 加这一行：重新分配后旧缓存失效
        }
        dno.ensure_size(N, DNO_M, DNO_Mb);
    }
};

// ---------------------------------------------------------------------
// 自适应步长积分主循环：从 t0 积分到 t_end
//
//   - 当剩余时间不足 dt 时，把 dt 收缩到恰好落在 t_end；
//   - 返回实际接受的步数；
//   - 若 dt 下溢到 阈值 以下或步数超过 max_steps，抛异常。
// 新增参数：
//   recorder  可选输出记录器；为 nullptr 或 disabled 时不做任何记录。
//             若提供，本函数会：
//               - 在 t0 处写出初始记录；
//               - 每一步接受后若已到达记录时刻，则写出记录。
//             调用者需保证 recorder 在调用前已用 (cfg, t0) 构造。
// ---------------------------------------------------------------------
size_t dno_time_integrate(
    DNOState& state,
    double t0, double t_end,
    double dt_init,
    double rtol, double atol,
    const DNOSystemParams& p,
    size_t max_steps = 100000,
    DNOOutputRecorder* recorder = nullptr);

#endif // DNO_TIME_INTEGRAL_H