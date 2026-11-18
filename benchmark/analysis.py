"""Analysis of the benchmark results (Paper B, Sec. VI).

Reads results/raw.csv (per seed) and results/contact_hist.csv (written by
collect_results.m) and computes
  * means and 95% confidence half-widths over seeds,
  * a contact-based predictor of the D2D offloading ratio that uses only the
    contact-length histogram of a model, and its error against simulation.

The predictor. For a pair of users, let h_l be the number of complete
contacts of length l seconds per pair over a period of T seconds. Treating
contacts as a stationary sequence, the probability that the pair offers at
least Dd consecutive seconds of contact inside a window of Td seconds that
starts at a random instant is approximately

    g(Dd) = [ sum_l h_l * ( max(l - Dd + 1, 0) + 1{l >= Dd} * (Td - Dd) ) ] / T

(first term: already in contact with enough residual time; second term: a
long-enough contact starts early enough in the window). With q_f the
probability that a given other user caches file f, and holders treated as
independent,

    offload = sum_f p_f [ q_f + (1 - q_f) * (1 - (1 - g q_f)^(U-1)) ]

(the requester holds f with probability q_f; the number of other holders is
Binomial(U-1, q_f)).
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
U, N, GAMMA, C_USER, TD, TSIM = 200, 1000, 0.8, 10, 60, 3600
T95 = {2: 12.706, 3: 4.303, 5: 2.776, 10: 2.262}


def zipf(n=N, g=GAMMA):
    w = np.arange(1, n + 1, dtype=float) ** (-g)
    return w / w.sum()


def inclusion_prob(p, c=C_USER, reps=20000, seed=1):
    """Probability that a user holds each file under sequential sampling of c
    distinct files with probability proportional to popularity (Monte Carlo,
    same procedure as run_one.m)."""
    rng = np.random.default_rng(seed)
    cnt = np.zeros(p.size)
    for _ in range(reps):
        q = p.copy()
        for _ in range(c):
            f = rng.choice(p.size, p=q / q.sum())
            cnt[f] += 1
            q[f] = 0.0
    return cnt / reps


def pair_prob(hist, n_pairs, dd, td=TD, T=TSIM):
    l = np.arange(1, hist.size + 1)
    resid = np.maximum(l - dd + 1, 0)
    start = (l >= dd) * max(td - dd, 0)
    return float((hist * (resid + start)).sum() / n_pairs / T)


def predict_offload(g, p, q, u=U):
    # the number of other holders is Binomial(u-1, q_f): E[(1-g)^H] = (1-g q_f)^(u-1)
    return float((p * (q + (1 - q) * (1 - (1 - g * q) ** (u - 1)))).sum())


def summarize(df, cols):
    gb = df.groupby(["model", "speed"])
    m = gb[cols].mean()
    n = gb.size().iloc[0]
    ci = gb[cols].std(ddof=1) * T95.get(n, 2.0) / np.sqrt(n)
    return m, ci, n




def load_all():
    raw = pd.read_csv(os.path.join(RES, "raw.csv"))
    h = pd.read_csv(os.path.join(RES, "contact_hist.csv"), header=None)
    h.columns = ["model", "speed", "seed", "pairs"] + list(range(1, h.shape[1] - 3))
    return raw, h


def d2d_prediction(raw, h, dds=(2, 5, 10, 20, 30)):
    """Predicted offloading per run, from that run's contact histogram."""
    p = zipf()
    q = inclusion_prob(p, reps=40000)
    rows = []
    for (_, r), (_, hr) in zip(raw.iterrows(), h.iterrows()):
        assert r.model == hr.model and r.speed == hr.speed and r.seed == hr.seed
        hist = hr.iloc[4:].to_numpy(dtype=float)
        for dd in dds:
            if r.model == "static":
                g = r.contact_frac          # pairs in range stay in range
            else:
                g = pair_prob(hist, hr.pairs, dd)
            rows.append(dict(model=r.model, speed=r.speed, seed=r.seed, Dd=dd, g=g,
                             pred=predict_offload(g, p, q), sim=r["d2d_Dd%d" % dd]))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    raw, h = load_all()
    D = d2d_prediction(raw, h)
    D.to_csv(os.path.join(RES, "d2d_prediction.csv"), index=False)
    M = D.groupby(["model", "speed", "Dd"])[["pred", "sim"]].mean()
    e = (M.pred - M.sim).abs()
    print("configurations: %d, mean |error| %.4f, 90%% below %.4f, max %.4f"
          % (len(e), e.mean(), e.quantile(0.9), e.max()))
