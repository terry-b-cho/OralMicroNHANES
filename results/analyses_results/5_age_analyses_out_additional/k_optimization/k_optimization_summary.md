# k optimization summary

Generated: 2026-06-08 03:57:08

Search range: k = 2..15; n_genera = 93; B_GAP = 100; B_STABILITY = 100.

**Chosen number of clusters: k = 9.**

## Evidence by panel

Each panel of `k_optimization_panels.pdf` shows the data behind a complementary criterion. The reading at k = 9 is summarized below.

| Panel | Criterion | Reading at k = 9 |
|---|---|---|
| A | Silhouette: local maxima | s_bar(9) = 0.253; k=9 is a local maximum |
| B | Gap statistic (Tibshirani SE-rule) | decision_signal(9) = 0.011; k=9 satisfies Gap(k) >= Gap(k+1) - SE(k+1) |
| C | W (within-cluster dispersion) and its derivatives | W''(9) = 0.645 (eps = 0.671); k=9 is the W''->0 plateau |
| D | Resolution ceiling (first-touch-the-floor) | min_size(9) = 3 (= phi = 3); k=9 first touches the resolution floor |
| E | Stability (bootstrap consensus) | mean_jaccard(9) = 0.899 (>= 0.75); k=9 is in the stable plateau |

Bootstrap consensus heatmaps for k = 2..15 are in `consensus_stability.pdf`; the k = 9 panel is boxed in teal.

