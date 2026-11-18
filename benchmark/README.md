# Paper B benchmark: user mobility models in edge caching

MATLAB code for Section VI of "Does User Mobility Help or Hurt Edge Caching?
A Survey and a Unified Benchmark". It runs unchanged in
GNU Octave (tested with Octave 8 on Ubuntu).

## Run

```matlab
run_sweep(1,1)        % 7 models x 6 speeds x 10 seeds -> results/main_*.mat
run_mismatch(10)      % train on one model, test on another (Table VII uses v = 10 m/s;
                      % run_mismatch(1) gives the walking-speed check)
run_mixed()           % 70% static users + 30% random walkers (Table VIII)
collect_results()     % summary.mat, summary.csv, raw.csv (per seed), contact_hist.csv
d2d_predict()         % contact-based prediction of D2D offloading, Eqs. (4)-(5)
run_rwp_warmup(5)     % classical vs stationary random waypoint (quoted in Sec. VI-A)
spatial_check()       % border and center occupancy of each model
plot_results()        % Figs. 2 and 3 (MATLAB)
```

`run_sweep(k, n)` splits the sweep across `n` processes, and
`run_mismatch(5, seeds)` plus `merge_mismatch(5)` does the same for the
mismatch test. The 10-seed sweep took about 70 minutes on two cores with Octave.

The figures in the paper were rendered with `plot_results.py` (matplotlib), which
reads `results/raw.csv` and `results/d2d_prediction.csv` (from `analysis.py`)
and draws the same panels as `plot_results.m`.

## Files

| File | Role |
|---|---|
| `bench_config.m` | All parameters (Table V of the paper) |
| `gen_mobility.m` | Seven mobility models; reflection at the walls; warm-up discarded |
| `sbs_assoc.m` | Nearest-SBS association on a 5 x 5 grid |
| `estimate_reach.m` | Reach matrix W used by the mobility-aware placement |
| `place_mobility_aware.m` | Greedy placement for objective (3) |
| `gen_requests.m`, `zipf_pmf.m` | Poisson requests, Zipf popularity |
| `eval_requests.m` | Snapshot and effective hit ratios, LRU, D2D offloading |
| `mobility_stats.m` | Speed, sojourn, handovers, reach, stay probability, contact statistics |
| `run_one.m` | One (model, speed, seed) point |
| `run_mixed.m` | Mixed population: share `P.mix_static` of users static |
| `d2d_predict.m`, `analysis.py` | Usable-contact probability g(Dd) and predicted offloading |

## Changes from the first draft

* Gauss-Markov users were clamped to the wall, which pinned fast users at the
  border; they are now reflected with a mirrored heading.
* Random-walk and Levy flights were drawn as straight lines between reflected
  endpoints, which cut corners at high speed; flights are now sampled every
  second along the reflected path.
* Seeds were raised from 3 to 10, and five D2D transfer times are evaluated.

To use a trace, replace the call to `gen_mobility` in `run_one.m` with a loader
that returns U-by-T position matrices sampled every second.
