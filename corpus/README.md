# Coded caching studies

`coded_studies.csv` is the coding sheet behind Table III of the paper. Each row is
one caching study, identified by its key in `paper/refs.bib`.

| Column | Values |
|---|---|
| `level` | how the caching decision uses mobility: parametric, prediction-driven, oblivious, not applicable |
| `mobility_swept` | Y: the evaluation varies a mobility parameter; N: it does not; ?: not yet checked; --: not applicable |
| `reported_effect` | effect of mobility on the gain as reported by the authors; n.c.: not coded |

Rows were coded from the full text where available and otherwise from the abstract.
The corpus is a draft and will be extended before submission.
