# Coded caching studies

`coded_studies.csv` is the coding sheet behind Table III of the paper. Each row is
one caching study, identified by its key in `paper/refs.bib`.

| Column | Values |
|---|---|
| `level` | how the caching decision uses mobility: parametric, prediction-driven, oblivious, not applicable |
| `mobility_swept` | Y: the evaluation varies a mobility parameter; N: it does not; A: only the abstract was available and it reports no sweep; --: not applicable |
| `reported_effect` | effect of mobility on the gain as reported by the authors; n.r.: not reported |

Rows were coded from the full text; the three rows marked A were coded from the abstract.
