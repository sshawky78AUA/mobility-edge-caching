"""Render the Section VI figures from the MATLAB benchmark output.

Inputs : results/raw.csv            (one row per model, speed and seed)
         results/d2d_prediction.csv (written by analysis.py)
Outputs: ../figures/fig_sbs_hit.pdf, ../figures/fig_d2d.pdf

Points are means over seeds; error bars are 95% confidence intervals
(Student t over seeds). The equivalent MATLAB plotting code is plot_results.m.
Palette: validated categorical slots 1-6; the static model is the gray reference.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
OUT = os.path.join(HERE, "..", "figures")
MODELS = ["static", "randwalk", "rwp", "gaussmarkov", "markovcell", "manhattan", "levy"]
NAMES = ["Static", "Random walk", "RWP", "Gauss-Markov", "Markov cell", "Manhattan", "Lévy walk"]
COLORS = ["#8a8a85", "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
MARKERS = ["", "o", "s", "^", "D", "v", "P"]
SPEEDS = [0.5, 1, 2, 5, 10, 20]
T95 = {3: 4.303, 5: 2.776, 10: 2.262}

plt.rcParams.update({
    "font.family": "serif", "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.edgecolor": "#6b6b66", "axes.linewidth": 0.6, "grid.color": "#e2e1dc",
    "grid.linewidth": 0.5, "pdf.fonttype": 42,
})


def stats(raw, key):
    g = raw.groupby(["model", "speed"])[key]
    n = raw.groupby(["model", "speed"]).size().iloc[0]
    return g.mean(), g.std(ddof=1) * T95.get(n, 2.0) / np.sqrt(n)


def clean(ax):
    ax.grid(True, which="major")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def speed_axis(ax):
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(SPEEDS))
    ax.xaxis.set_major_formatter(FixedFormatter([("%g" % v) for v in SPEEDS]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlim(0.4, 25)
    ax.set_xlabel("Travel speed $v$ (m/s)")
    clean(ax)


def draw_models(ax, raw, key):
    mean, ci = stats(raw, key)
    handles = []
    for m, n, c, mk in zip(MODELS, NAMES, COLORS, MARKERS):
        y = np.array([mean[(m, v)] for v in SPEEDS])
        e = np.array([ci[(m, v)] for v in SPEEDS])
        if m == "static":
            h, = ax.plot(SPEEDS, y, color=c, lw=1.2, ls="--", label=n)
        else:
            h = ax.errorbar(SPEEDS, y, yerr=e, color=c, lw=1.4, marker=mk, ms=4.2,
                            markeredgecolor="white", markeredgewidth=0.6,
                            elinewidth=0.8, capsize=1.6, label=n)
        handles.append(h)
    return handles


def legend(fig, handles):
    fig.legend(handles, NAMES, loc="outside lower center", ncol=7, frameon=False,
               handlelength=2.2, columnspacing=1.2)


def fig_sbs(raw):
    fig, axs = plt.subplots(2, 2, figsize=(7.16, 4.0), constrained_layout=True)
    axs = axs.ravel()
    mpc = raw.eff_mpc.mean()
    panels = [("eff_ma", "(a) Mobility-aware placement, effective hit ratio", (0.32, 0.54)),
              ("snap_ma", "(b) Mobility-aware placement, snapshot hit ratio", (0.32, 0.54)),
              ("eff_lru", "(c) LRU at the serving SBS, effective hit ratio", (0.04, 0.26))]
    for ax, (key, title, yl) in zip(axs[:3], panels):
        handles = draw_models(ax, raw, key)
        if key != "eff_lru":
            ax.axhline(mpc, color="#2b2b28", lw=0.8, ls=":")
            ax.text(0.45, mpc + 0.004, "Most popular", fontsize=6.5, color="#2b2b28")
        ax.set_title(title, loc="left"); ax.set_ylim(*yl)
        ax.set_ylabel("Hit ratio")
        speed_axis(ax)
    # (d) gain against the measured reach
    ax = axs[3]
    m = raw.groupby(["model", "speed"]).mean(numeric_only=True).reset_index()
    m["gain"] = 100 * (m.eff_ma / m.eff_mpc - 1)
    for mod, c, mk in zip(MODELS, COLORS, MARKERS):
        d = m[m.model == mod]
        ax.plot(d.mean_reach, d.gain, ls="none", marker=mk if mk else "o", ms=4.5,
                color=c, markeredgecolor="white", markeredgewidth=0.6)
    mv = m[(m.model != "static") & (m.model != "rwp")]
    k = (mv.gain * (mv.mean_reach - 1)).sum() / ((mv.mean_reach - 1) ** 2).sum()
    xs = np.linspace(1, 3.4, 50)
    ax.plot(xs, k * (xs - 1), color="#2b2b28", lw=0.8, ls=":")
    ax.text(1.05, 23.5, "Dotted: %.1f%% per additional SBS" % k, fontsize=6.5, color="#2b2b28")
    ax.set_xlim(0.95, 3.4); ax.set_ylim(-1, 27)
    ax.set_xlabel("Reach: SBSs usable within the deadline (mean)")
    ax.set_ylabel("Gain (%)")
    ax.set_title("(d) Gain over most popular against reach", loc="left")
    clean(ax)
    legend(fig, handles)
    fig.savefig(os.path.join(OUT, "fig_sbs_hit.pdf"))
    plt.close(fig)
    return k


def fig_d2d(raw, P):
    fig, axs = plt.subplots(2, 2, figsize=(7.16, 4.0), constrained_layout=True)
    axs = axs.ravel()
    for ax, dd, lab, yl in zip(axs[:3], [2, 10, 30], "abc",
                               [(0.14, 0.50), (0.11, 0.27), (0.07, 0.19)]):
        handles = draw_models(ax, raw, "d2d_Dd%d" % dd)
        ax.set_title("(%s) Transfer time $D_d$ = %d s" % (lab, dd), loc="left")
        ax.set_ylabel("D2D offloading ratio"); ax.set_ylim(*yl)
        speed_axis(ax)
    ax = axs[3]
    M = P.groupby(["model", "speed", "Dd"])[["g", "sim", "pred"]].mean().reset_index()
    for mod, c, mk in zip(MODELS, COLORS, MARKERS):
        d = M[M.model == mod]
        ax.plot(d.g, d.sim, ls="none", marker=mk if mk else "o", ms=3.8, color=c,
                markeredgecolor="white", markeredgewidth=0.5)
    gs = np.logspace(-3.3, -0.7, 80)
    import sys
    sys.path.insert(0, HERE)
    from analysis import zipf, inclusion_prob, predict_offload
    p = zipf(); q = inclusion_prob(p, reps=40000)
    ax.plot(gs, [predict_offload(g, p, q) for g in gs], color="#2b2b28", lw=1.0)
    ax.set_xscale("log"); ax.set_xlim(5e-4, 0.2); ax.set_ylim(0.07, 0.52)
    ax.set_xlabel(r"Usable-contact probability per pair, $g(D_d)$")
    ax.set_ylabel("D2D offloading ratio")
    ax.set_title("(d) All runs against $g(D_d)$; solid line: prediction", loc="left")
    clean(ax)
    legend(fig, handles)
    fig.savefig(os.path.join(OUT, "fig_d2d.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    raw = pd.read_csv(os.path.join(RES, "raw.csv"))
    P = pd.read_csv(os.path.join(RES, "d2d_prediction.csv"))
    k = fig_sbs(raw)
    print("gain per extra reachable SBS (%%, excluding RWP): %.2f" % k)
    fig_d2d(raw, P)
