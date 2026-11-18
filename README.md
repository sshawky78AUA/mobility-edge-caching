# User Mobility in Edge Caching

Code, data and LaTeX source for

> S. Hosny and J. Tadrous, "User Mobility in Edge Caching: A Survey of Models,
> Their Effect on Caching Gain, and a Unified Benchmark," manuscript in preparation for
> *IEEE Communications Surveys & Tutorials*, 2026.

## Contents

| Folder | What it holds |
|---|---|
| `benchmark/` | MATLAB benchmark of Section VI (runs in MATLAB and GNU Octave), Python plotting and analysis, and every result file |
| `benchmark/results/` | Per-run results (`main_*.mat`), per-seed table (`raw.csv`), contact-length histograms, mismatch and warm-up runs |
| `corpus/` | Coding sheet of the 24 caching studies (Table III) |
| `paper/` | LaTeX source (IEEEtran journal class), bibliography and figures |

## Reproducing the results

In MATLAB or Octave, from `benchmark/`:

```matlab
run_sweep(1,1)      % 7 models x 6 speeds x 10 seeds (about 70 min on two cores in Octave)
run_mismatch(10)    % train/test mismatch at 10 m/s (and run_mismatch(1) at walking speed)
run_mixed()         % mixed population of static users and random walkers
collect_results()   % summary.csv, raw.csv, contact_hist.csv
d2d_predict()       % contact-based prediction of D2D offloading, Eqs. (4)-(5)
plot_results()      % Figs. 2 and 3
```

The paper figures were rendered with `python3 analysis.py` followed by
`python3 plot_results.py` (numpy, pandas, scipy, matplotlib), which read the same
result files. See `benchmark/README.md` for the file-by-file description.

## Building the paper

```
cd paper
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## License

Code: MIT (`LICENSE`). Results and coding sheet: CC BY 4.0 (`LICENSE-DATA.md`).
