"""Render the Section VI figures from the MATLAB benchmark output.

Inputs : results/raw.csv            (one row per model, speed and seed)
         results/contact_hist.csv   (contact-length histograms, collect_results.m)
         results/d2d_prediction.csv (written by analysis.py)
Outputs: ../paper/figures/fig_sbs_hit.pdf, fig_d2d.pdf, fig_predict.pdf

Each figure answers one question per panel. Curves are means over seeds and the
shaded bands are 95% confidence intervals (Student t over seeds). The
equivalent MATLAB plotting code is plot_results.m.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analysis import zipf, inclusion_prob, predict_offload  # noqa: E402

RES = os.path.join(HERE, "results")
OUT = next(d for d in (os.path.join(HERE, "..", "paper", "figures"),
                      os.path.join(HERE, "..", "figures")) if os.path.isdir(d))
MODELS = ["static", "randwalk", "rwp", "gaussmarkov", "markovcell", "manhattan", "levy"]
NAMES = ["Static", "Random walk", "RWP", "Gauss-Markov", "Markov cell", "Manhattan", "Lévy walk"]
COLORS = ["#8a8a85", "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
MARKERS = ["o", "o", "s", "^", "D", "v", "P"]
DDS = [2, 5, 10, 20, 30]
DD_COLORS = ["#9ec5f0", "#5e9fe6", "#2a78d6", "#1a55a3", "#0d3266"]   # one hue, light to dark
TICKS = [0.5, 1, 2, 5, 10, 20]
T95 = {2: 12.706, 3: 4.303, 5: 2.776, 10: 2.262}
INK, MUTED = "#2b2b28", "#6b6b66"

plt.rcParams.update({
    "font.family": "serif", "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "grid.color": "#e2e1dc",
    "grid.linewidth": 0.5, "pdf.fonttype": 42,
})


def ci_stats(raw, key):
    g = raw.groupby(["model", "speed"])[key]
    n = raw.groupby(["model", "speed"]).size().iloc[0]
    return g.mean(), g.std(ddof=1) * T95.get(n, 2.0) / np.sqrt(n)


def clean(ax):
    ax.grid(True, which="major")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def speed_axis(ax, label=True):
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(TICKS))
    ax.xaxis.set_major_formatter(FixedFormatter([("%g" % v) for v in TICKS]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlim(0.45, 22)
    if label:
        ax.set_xlabel("Travel speed $v$ (m/s)")
    clean(ax)


def speed_regimes(ax, y):
    """Light labels for walking, cycling and driving speeds."""
    for lo, hi, txt in [(0.45, 2, "walking"), (3, 7, "cycling"), (10, 22, "driving")]:
        ax.axvspan(lo, hi, color="#f3f2ee", zorder=0, lw=0)
        ax.text(np.sqrt(lo * hi), y, txt, ha="center", va="top", fontsize=6.3, color=MUTED)


def band(ax, x, m, e, color, label=None, ls="-", lw=1.4):
    ax.fill_between(x, m - e, m + e, color=color, alpha=0.18, lw=0)
    return ax.plot(x, m, color=color, lw=lw, ls=ls, label=label)[0]


def model_means(raw):
    m = raw.groupby(["model", "speed"]).mean(numeric_only=True).reset_index()
    m["gain_eff"] = 100 * (m.eff_ma / m.eff_mpc - 1)
    m["gain_snap"] = 100 * (m.snap_ma / m.snap_mpc - 1)
    return m


# ---------------------------------------------------------------- Fig. 2 (SBS)
def fig_sbs(raw):
    fig, axs = plt.subplots(1, 3, figsize=(7.16, 2.55), constrained_layout=True)
    axs = axs.ravel()
    raw = raw.copy()
    raw["gain_eff"] = 100 * (raw.eff_ma / raw.eff_mpc - 1)
    mean, ci = ci_stats(raw, "gain_eff")
    speeds = sorted(raw.speed.unique())
    M = model_means(raw)

    # (a) gain against speed
    ax = axs[0]
    speed_regimes(ax, 27.5)
    handles = []
    for mod, nm, c in zip(MODELS, NAMES, COLORS):
        y = np.array([mean[(mod, v)] for v in speeds])
        e = np.array([ci[(mod, v)] for v in speeds])
        handles.append(band(ax, speeds, y, e, c, nm, ls="--" if mod == "static" else "-",
                            lw=1.1 if mod == "static" else 1.4))
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_ylim(-2, 28)
    ax.set_ylabel("Gain over most popular (%)")
    ax.set_title("(a) Gain over MPC against speed", loc="left")
    speed_axis(ax)

    # (b) snapshot gain against effective gain
    ax = axs[1]
    mv = M[M.model != "static"]
    ax.axvspan(-30, 0, ymin=0, ymax=1, color="#fbeee8", lw=0, zorder=0)
    ax.text(-25, 3.5, "snapshot: loss\ndelivery: gain", fontsize=6.5,
            color="#b04a1f", va="bottom")
    for mod, c, mk in zip(MODELS[1:], COLORS[1:], MARKERS[1:]):
        d = mv[mv.model == mod]
        ax.plot(d.gain_snap, d.gain_eff, ls="none", marker=mk, ms=4.2, color=c,
                markeredgecolor="white", markeredgewidth=0.5)
    ax.plot([-30, 30], [-30, 30], color=MUTED, lw=0.7, ls=":")
    ax.axhline(0, color=INK, lw=0.6); ax.axvline(0, color=INK, lw=0.6)
    ax.set_xlim(-26, 3); ax.set_ylim(-2, 28)
    ax.set_xlabel("Gain, snapshot hit ratio (%)")
    ax.set_ylabel("Gain, effective hit ratio (%)")
    ax.set_title("(b) Snapshot vs effective gain", loc="left")
    clean(ax)

    # (c) LRU: effective hit ratio against the probability of finishing in the first cell
    ax = axs[2]
    snap = raw.snap_lru.mean()
    for mod, c, mk in zip(MODELS, COLORS, MARKERS):
        d = M[M.model == mod]
        ax.plot(d.p_stay_first, d.eff_lru, ls="none", marker=mk, ms=4.2, color=c,
                markeredgecolor="white", markeredgewidth=0.5)
    xs = np.linspace(0, 1, 20)
    ax.plot(xs, snap * xs, color=INK, lw=0.8)
    ax.text(0.98, 0.06, "line: snapshot ratio\n" + r"$\times$ P(stay $\geq D_s$)",
            fontsize=6.5, color=INK, ha="right")
    ax.set_xlim(0.2, 1.02); ax.set_ylim(0.04, 0.26)
    ax.set_xlabel(r"P(stay $\geq D_s$ in the first cell)")
    ax.set_ylabel("LRU effective hit ratio")
    ax.set_title("(c) LRU against handover risk", loc="left")
    clean(ax)

    hs = [plt.Line2D([], [], color=c, lw=1.4, marker=mk, ms=4, ls="--" if i == 0 else "-",
                     markeredgecolor="white", markeredgewidth=0.5)
          for i, (c, mk) in enumerate(zip(COLORS, MARKERS))]
    fig.legend(hs, NAMES, loc="outside lower center", ncol=7, frameon=False,
               handlelength=2.4, columnspacing=1.2)
    fig.savefig(os.path.join(OUT, "fig_sbs_hit.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- Fig. 3 (D2D)
def fig_d2d(raw, hist):
    fig, axs = plt.subplots(1, 3, figsize=(7.16, 2.55), constrained_layout=True)
    axs = axs.ravel()
    speeds = sorted(raw.speed.unique())
    base = raw[raw.model == "static"]
    raw = raw.copy()

    # (a) relative change against speed, random walk, every transfer time
    ax = axs[0]
    speed_regimes(ax, 215)
    rw = raw[raw.model == "randwalk"]
    hs = []
    for dd, c in zip(DDS, DD_COLORS):
        key = "d2d_Dd%d" % dd
        s0 = base[key].mean()
        rel = rw.assign(rel=100 * (rw[key] / s0 - 1))
        m, e = ci_stats(rel, "rel")
        y = np.array([m[("randwalk", v)] for v in speeds])
        ee = np.array([e[("randwalk", v)] for v in speeds])
        hs.append(band(ax, speeds, y, ee, c, "$D_d$ = %d s" % dd))
    ax.axhline(0, color=INK, lw=0.6)
    ax.text(21, 4, "static network", fontsize=6.3, color=INK, ha="right")
    ax.set_ylim(-60, 220)
    ax.set_ylabel("Offloading vs static (%)")
    ax.set_title("(a) Offloading vs static, random walk", loc="left")
    ax.legend(handles=hs, loc="upper left", bbox_to_anchor=(0.0, 0.9), frameon=False,
              handlelength=1.6, labelspacing=0.25)
    speed_axis(ax)

    # (b) mechanism: share of contacts that last at least D_d
    ax = axs[1]
    h = hist[hist.model == "randwalk"]
    L = np.arange(1, h.shape[1] - 3)
    for v, c, ls in zip([1, 5, 20], ["#9ec5f0", "#2a78d6", "#0d3266"], ["-", "-", "-"]):
        cnt = h[h.speed == v].iloc[:, 4:].to_numpy(float).sum(axis=0)
        ccdf = cnt[::-1].cumsum()[::-1] / cnt.sum()
        ax.plot(L, ccdf, color=c, lw=1.4, ls=ls, label="$v$ = %g m/s" % v)
    for dd in DDS:
        ax.axvline(dd, color=MUTED, lw=0.5, ls=":")
        ax.text(dd, 1.12, "%d" % dd, fontsize=6.3, color=MUTED,
                ha={20: "right", 30: "left"}.get(dd, "center"))
    ax.text(60, 1.12, "$D_d$ (s)", ha="left", fontsize=6.3, color=MUTED)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1, 600); ax.set_ylim(3e-3, 1.0)
    ax.set_xlabel("Contact duration $l$ (s)")
    ax.set_ylabel(r"Share of contacts lasting $\geq l$")
    ax.set_title("(b) Contact durations, random walk", loc="left", pad=11)
    ax.legend(loc="lower left", frameon=False, handlelength=1.6, labelspacing=0.25)
    clean(ax)

    # (c) the best speed for each transfer time and model (Levy walk excluded:
    #     its offloading barely changes with speed, so the maximum is ill defined)
    ax = axs[2]
    B = best_speeds(raw)
    for mod, c, mk in zip(MODELS[1:6], COLORS[1:6], MARKERS[1:6]):
        d = B[B.model == mod]
        ax.plot(d.Dd, d.best_speed, color=c, lw=1.2)
        for _, r in d.iterrows():
            ax.plot(r.Dd, r.best_speed, marker=mk, ms=4.2, color=c, ls="none",
                    markerfacecolor=c if r.interior else "white",
                    markeredgecolor="white" if r.interior else c, markeredgewidth=0.6)
    fit = B[(B.model != "levy") & B.interior]
    kap = np.exp(np.mean(np.log(fit.best_speed * fit.Dd)))
    xs = np.linspace(1.8, 34, 50)
    ax.plot(xs, kap / xs, color=INK, lw=0.8, ls=":")
    ax.text(2.1, 0.55, r"dotted: $v^\star D_d$ = %.0f m" % kap, fontsize=6.5, color=INK)
    ax.text(2.1, 0.9, "hollow: optimum at the\nedge of the speed range", fontsize=6.3,
            color=MUTED, va="bottom")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.xaxis.set_major_locator(FixedLocator(DDS))
    ax.xaxis.set_major_formatter(FixedFormatter([str(d) for d in DDS]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_locator(FixedLocator(TICKS))
    ax.yaxis.set_major_formatter(FixedFormatter([("%g" % v) for v in TICKS]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xlim(1.7, 36); ax.set_ylim(0.4, 26)
    ax.set_xlabel("Transfer time $D_d$ (s)")
    ax.set_ylabel(r"Best speed $v^\star$ (m/s)")
    ax.set_title("(c) Best speed per transfer time", loc="left")
    clean(ax)

    hs = [plt.Line2D([], [], color=c, lw=1.2, marker=mk, ms=4,
                     markeredgecolor="white", markeredgewidth=0.5)
          for c, mk in zip(COLORS[1:6], MARKERS[1:6])]
    fig.legend(hs, NAMES[1:6], loc="outside lower center", ncol=5,
               frameon=False, handlelength=2.4, columnspacing=1.2)
    fig.savefig(os.path.join(OUT, "fig_d2d.pdf"))
    plt.close(fig)


# ------------------------------------------- Fig. 4 (what predicts the outcome)
def fig_predict(raw, P):
    fig, axs = plt.subplots(1, 2, figsize=(7.16, 2.45), constrained_layout=True)
    M = model_means(raw)
    # (a) SBS: gain against the measured reach
    ax = axs[0]
    for mod, c, mk in zip(MODELS, COLORS, MARKERS):
        d = M[M.model == mod]
        ax.plot(d.mean_reach, d.gain_eff, ls="none", marker=mk, ms=4.2, color=c,
                markeredgecolor="white", markeredgewidth=0.5)
    fit = M[(M.model != "static") & (M.model != "rwp")]
    k = (fit.gain_eff * (fit.mean_reach - 1)).sum() / ((fit.mean_reach - 1) ** 2).sum()
    r = np.corrcoef(fit.mean_reach, fit.gain_eff)[0, 1]
    xs = np.linspace(1, 3.4, 50)
    ax.plot(xs, k * (xs - 1), color=INK, lw=0.8, ls=":")
    ax.text(1.03, 26, "%.1f%% per additional SBS ($r$ = %.2f)" % (k, r), fontsize=6.5,
            color=INK, va="top")
    ax.annotate("RWP: reach spent\nin central cells", xy=(2.1, 7.0), xytext=(2.45, 1.5),
                fontsize=6.3, color="#b04a1f",
                arrowprops=dict(arrowstyle="-", color="#b04a1f", lw=0.6))
    ax.set_xlim(0.95, 3.4); ax.set_ylim(-2, 28)
    ax.set_xlabel("Reach: SBSs usable within the deadline (mean)")
    ax.set_ylabel("Gain over most popular (%)")
    ax.set_title("(a) SBS caching: reach predicts the gain", loc="left")
    clean(ax)

    # (b) D2D: prediction against simulation for every configuration
    ax = axs[1]
    Q = P.groupby(["model", "speed", "Dd"])[["sim", "pred"]].mean().reset_index()
    for dd, c in zip(DDS, DD_COLORS):
        d = Q[Q.Dd == dd]
        ax.plot(d.pred, d.sim, ls="none", marker="o", ms=3.4, color=c,
                markeredgecolor="white", markeredgewidth=0.4, label="$D_d$ = %d s" % dd)
    ax.plot([0, 0.6], [0, 0.6], color=INK, lw=0.8)
    err = (Q.pred - Q.sim).abs()
    ax.text(0.075, 0.47, "%d configurations\nmean |error| %.3f, max %.3f"
            % (len(Q), err.mean(), err.max()), fontsize=6.5, color=INK, va="top")
    ax.set_xlim(0.07, 0.52); ax.set_ylim(0.07, 0.52)
    ax.set_xlabel("Predicted offloading ratio, Eq. (5)")
    ax.set_ylabel("Simulated offloading ratio")
    ax.set_title("(b) D2D caching: contact statistics predict offloading", loc="left")
    ax.legend(loc="lower right", frameon=False, handlelength=1.0, labelspacing=0.2)
    clean(ax)

    hs = [plt.Line2D([], [], color=c, lw=0, marker=mk, ms=4.5,
                     markeredgecolor="white", markeredgewidth=0.5)
          for c, mk in zip(COLORS, MARKERS)]
    axs[0].legend(hs, NAMES, loc="upper left", bbox_to_anchor=(0.0, 0.9), frameon=False,
                  handletextpad=0.2,
                  labelspacing=0.2, fontsize=6.5)
    fig.savefig(os.path.join(OUT, "fig_predict.pdf"))
    plt.close(fig)
    return k, r


def best_speeds(raw):
    """Speed that maximizes the offloading ratio, refined by a parabola in log v
    through the best grid point and its two neighbours."""
    base = raw[raw.model == "static"]
    M = raw.groupby(["model", "speed"]).mean(numeric_only=True).reset_index()
    rows = []
    for mod in MODELS[1:]:
        d = M[M.model == mod].sort_values("speed")
        x = np.log(d.speed.to_numpy())
        for dd in DDS:
            y = d["d2d_Dd%d" % dd].to_numpy()
            i = int(np.argmax(y))
            interior = 0 < i < len(y) - 1
            xb = x[i]
            if interior:
                a, b, _ = np.polyfit(x[i - 1:i + 2], y[i - 1:i + 2], 2)
                if a < 0:
                    xb = float(np.clip(-b / (2 * a), x[i - 1], x[i + 1]))
            rows.append((mod, dd, float(np.exp(xb)), y[i], base["d2d_Dd%d" % dd].mean(),
                         interior))
    return pd.DataFrame(rows, columns=["model", "Dd", "best_speed", "best", "static",
                                       "interior"])


if __name__ == "__main__":
    raw = pd.read_csv(os.path.join(RES, "raw.csv"))
    hist = pd.read_csv(os.path.join(RES, "contact_hist.csv"), header=None)
    hist.columns = ["model", "speed", "seed", "pairs"] + list(range(1, hist.shape[1] - 3))
    P = pd.read_csv(os.path.join(RES, "d2d_prediction.csv"))
    fig_sbs(raw)
    fig_d2d(raw, hist)
    k, r = fig_predict(raw, P)
    print("gain per extra reachable SBS (%%, excluding RWP): %.2f, r = %.3f" % (k, r))
    print(best_speeds(raw).to_string())
