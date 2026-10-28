#!/usr/bin/env python3
"""
dno_animate.py — 从 DNO 输出文件生成 η–x 与 β–x 的演化视频/动图（单 Y 轴）
                 并可选地在图例中显示每时刻的能量/动量（EVI.txt）

输入：
    eta.txt   自由表面高程
    beta.txt  变底部高程
    EVI.txt   能量/动量时间历程（可选）：每行 t,Ek,Ep,I（逗号分隔，无表头）
文件格式（eta.txt / beta.txt）：
    第 1 行：0, x0, x1, ..., x_{N-1}
    后续行：t_i, f(x0), f(x1), ..., f(x_{N-1})

显示变换：
    η_display = η * eta_scale
    β_display = (β - beta_offset) * beta_scale

默认：
    eta_scale   = 1.0
    beta_offset = β.min()   （把静态基底平移到 0）
    beta_scale  = 1.0
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# ---------------------------------------------------------------------
# 数据读取
# ---------------------------------------------------------------------
def load_wide_table(path: Path):
    """
    读取「宽表」格式：
        第 1 行  : 0, x0, x1, ..., x_{N-1}
        第 k 行  : t_k, f(x0), f(x1), ..., f(x_{N-1})
    返回 (x, t, F)，F.shape == (nt, N)。
    """
    data = np.loadtxt(path, delimiter=",")
    if data.ndim != 2 or data.shape[0] < 2:
        raise ValueError(f"{path}: 至少需要 1 行表头 + 1 行数据")
    x = data[0, 1:]
    t = data[1:, 0]
    F = data[1:, 1:]
    return x, t, F


def load_EVI(path: Path):
    """
    读取 EVI.txt：每行 t,Ek,Ep,I,V（5 列，逗号分隔，无表头）。
    返回 (t, Ek, Ep, I, V)，五者长度相同。
    """
    data = np.loadtxt(path, delimiter=",")
    if data.ndim != 2 or data.shape[1] != 5:
        raise ValueError(
            f"{path}: 期望 5 列 (t,Ek,Ep,I,V)，实际 shape={data.shape}")
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]


def padded_range(lo: float, hi: float, pad: float = 0.08):
    span = hi - lo
    if span == 0.0:
        span = max(abs(hi), 1.0)
    return lo - pad * span, hi + pad * span


# ---------------------------------------------------------------------
# 动画主函数
# ---------------------------------------------------------------------
def build_animation(
    eta_path: Path,
    beta_path: Path,
    out_path: Path,
    fps: int = 20,
    dpi: int = 150,
    eta_scale: float = 1.0,
    beta_offset: float | None = None,
    beta_scale: float = 1.0,
    evi_path: Path | None = None,
):
    # ---- 读 eta / beta ----
    x, t, eta = load_wide_table(eta_path)
    x_b, t_b, beta = load_wide_table(beta_path)

    if eta.shape != beta.shape:
        raise ValueError(
            f"eta.txt 与 beta.txt 形状不一致：{eta.shape} vs {beta.shape}"
        )
    if not np.allclose(t, t_b):
        raise ValueError("eta.txt 与 beta.txt 时间列不一致")

    nt = len(t)
    print(f"读取完成：{nt} 个时间步, {len(x)} 个空间点")
    print(f"  t ∈ [{t[0]:.6f}, {t[-1]:.6f}]")

    # ---- 读 EVI（可选）----
    E_total = None
    Ek = Ep = I_mom = V_mom = None
    if evi_path is not None:
        t_evi, Ek, Ep, I_mom, V_mom = load_EVI(evi_path)

        if len(t_evi) != len(t):
            warnings.warn(
                f"EVI.txt 行数({len(t_evi)}) 与 eta.txt 行数({len(t)}) 不一致，"
                "将按时间线性插值到 eta 的时间轴"
            )
            Ek    = np.interp(t, t_evi, Ek)
            Ep    = np.interp(t, t_evi, Ep)
            I_mom = np.interp(t, t_evi, I_mom)
            V_mom = np.interp(t, t_evi, V_mom)
        elif not np.allclose(t, t_evi, rtol=0, atol=1e-12):
            warnings.warn("EVI.txt 时间列与 eta.txt 时间列存在浮点差异")

        E_total = Ek + Ep
        print(f"  已读 EVI.txt："
              f"Ek ∈ [{Ek.min():.6e}, {Ek.max():.6e}], "
              f"Ep ∈ [{Ep.min():.6e}, {Ep.max():.6e}], "
              f"E  ∈ [{E_total.min():.6e}, {E_total.max():.6e}], "
              f"I  ∈ [{I_mom.min():.6e}, {I_mom.max():.6e}], "
              f"V  ∈ [{V_mom.min():.6e}, {V_mom.max():.6e}]")

    # ---- β 的显示变换 ----
    if beta_offset is None:
        beta_offset = float(beta.min())
    beta_disp = (beta - beta_offset) * beta_scale

    # 用于图例/坐标轴显示的 β 基线标签（处理符号，避免出现 "β − -1"）
    if beta_offset < 0:
        beta_label = f"(β + {abs(beta_offset):.6g}) * {beta_scale:.6g}"
    elif beta_offset > 0:
        beta_label = f"(β − {beta_offset:.6g}) * {beta_scale:.6g}"
    else:
        beta_label = f"β * {beta_scale:.6g}"

    # ---- η 的显示变换 ----
    eta_disp = eta * eta_scale
    eta_label = f"η * {eta_scale:.6g}"

    print(f"  η 显示：eta_scale = {eta_scale:.6g}")
    print(f"    η_display 范围 = [{eta_disp.min():.6e}, {eta_disp.max():.6e}]")
    print(f"  β 显示：beta_offset = {beta_offset:.6g}, "
          f"beta_scale = {beta_scale:.6g}")
    print(f"    β_display 范围 = [{beta_disp.min():.6e}, {beta_disp.max():.6e}]")

    # ---- 统一的 Y 轴范围 ----
    y_lo, y_hi = padded_range(
        float(min(eta_disp.min(), beta_disp.min())),
        float(max(eta_disp.max(), beta_disp.max())),
    )

    # ---- 画布与单 Y 轴 ----
    fig, ax = plt.subplots(figsize=(10.5, 5.8), dpi=dpi)

    c_eta  = "#1f77b4"     # 蓝
    c_beta = "#d62728"     # 红

    (line_eta,)  = ax.plot([], [], color=c_eta,  lw=2.3, label=eta_label)
    (line_beta,) = ax.plot([], [], color=c_beta, lw=2.3, label=beta_label)

    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y_lo, y_hi)

    ax.set_xlabel("x")
    if eta_scale != 1.0:
        ax.set_ylabel(f"η × {eta_scale:.3e}   /   {beta_label}")
    else:
        ax.set_ylabel(f"η   /   {beta_label}")
    ax.grid(True, ls=":", lw=0.6, alpha=0.6)
    ax.set_axisbelow(True)

    fig.tight_layout()

    # ---- 图例标题：t（+可选能量/动量），随帧更新 ----
    def format_title(k: int) -> str:
        s = f"t = {t[k]:.4f}"
        if E_total is not None:
            s += (f"\nE  = {E_total[k]:.8e}"
                  f"\nEk = {Ek[k]:.8e}"
                  f"\nEp = {Ep[k]:.8e}"
                  f"\nI  = {I_mom[k]:.8e}"
                  f"\nV  = {V_mom[k]:.8e}")
        return s

    legend_title0 = format_title(0)

    legend = ax.legend(
        [line_eta, line_beta], [eta_label, beta_label],
        loc="upper right",
        framealpha=0.9,
        title=legend_title0,
        title_fontsize=9,
    )

    # ---- 帧回调 ----
    def init():
        line_eta.set_data([], [])
        line_beta.set_data([], [])
        legend.set_title(legend_title0)
        ax.set_title("Wave / Bathymetry evolution")
        return line_eta, line_beta, legend

    def update(k: int):
        line_eta.set_data(x, eta_disp[k])
        line_beta.set_data(x, beta_disp[k])
        legend.set_title(format_title(k))
        ax.set_title(f"Wave / Bathymetry evolution   ({k+1}/{nt})")
        return line_eta, line_beta, legend

    ani = animation.FuncAnimation(
        fig, update,
        frames=nt,
        init_func=init,
        interval=1000.0 / fps,
        blit=False,
    )

    # ---- 保存 ----
    suffix = out_path.suffix.lower()
    if suffix == ".gif":
        ani.save(str(out_path), writer="pillow", fps=fps)
    elif suffix in (".mp4", ".mov", ".mkv", ".webm"):
        ani.save(str(out_path), writer="ffmpeg", fps=fps, dpi=dpi)
    else:
        raise ValueError(f"不支持的输出格式：{suffix}（请用 .mp4 或 .gif）")

    plt.close(fig)
    print(f"[OK] 已生成：{out_path}")


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="根据 DNO 输出文件生成 η–x 与 β–x 的演化动画（单 Y 轴，"
                    "可选叠加 EVI 能量/动量显示）"
    )
    ap.add_argument("--eta",  default="eta.txt",
                    help="η 数据文件（默认 eta.txt）")
    ap.add_argument("--beta", default="beta.txt",
                    help="β 数据文件（默认 beta.txt）")
    ap.add_argument("--evi",  default="EVI.txt",
                    help="EVI 数据文件（默认 EVI.txt；传 '' 表示禁用）")
    ap.add_argument("--out",  default="eta_beta_animation.mp4",
                    help="输出文件名；后缀决定格式（.mp4 / .gif）")
    ap.add_argument("--fps",  type=int, default=20,
                    help="帧率，默认 20")
    ap.add_argument("--dpi",  type=int, default=150,
                    help="渲染 dpi，默认 150")
    ap.add_argument("--eta-scale", type=float, default=1.0,
                    help="η 的显示缩放因子（仅用于可视化）。默认 1.0。")
    ap.add_argument("--beta-offset", type=float, default=None,
                    help="β 的显示基线（减法项）。默认取 β 的最小值。")
    ap.add_argument("--beta-scale", type=float, default=1.0,
                    help="β 的显示缩放因子。默认 1.0。")
    args = ap.parse_args()

    # 解析 EVI 文件（可选）
    evi: Path | None = None
    if args.evi:                      # 空字符串 → 禁用
        p = Path(args.evi)
        if p.exists():
            evi = p
        else:
            print(f"[warn] 未找到 {p}，将不显示能量信息")
            evi = None

    build_animation(
        Path(args.eta), Path(args.beta), Path(args.out),
        fps=args.fps, dpi=args.dpi,
        eta_scale=args.eta_scale,
        beta_offset=args.beta_offset,
        beta_scale=args.beta_scale,
        evi_path=evi,
    )


if __name__ == "__main__":
    main()