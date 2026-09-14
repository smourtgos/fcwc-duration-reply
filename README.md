# Replication package

**Interrogation Duration and the Estimation of False Confession Wrongful Conviction Risk: A Reply to Smith and Colleagues**
Scott M. Mourtgos and Ian T. Adams

## What this is

One script, `reply_analysis.R`, reproduces every number and every figure in the reply, in the order they appear in the text. It re-runs the interrogation-duration model of Smith et al. (2026) with the study inputs and model functions from their deposited code, changes only which studies enter the two duration likelihoods, and holds the FCWC base rate at their posterior median (0.019). No part of their model is re-specified.

Smith, T. B., Catlin, M., May, B., Redlich, A. D., Meissner, C. A., & Kelly, C. E. (2026). Simulating the probability of false confession-wrongful conviction conditional on features of police interrogation: A response to Mourtgos and Adams (2026). *Journal of Criminal Justice, 107*, 102742. https://doi.org/10.1016/j.jcrimjus.2026.102742

## Contents

| Path | What it is |
|---|---|
| `reply_analysis.R` | The analysis, start to finish. |
| `smith_et_al_2026_deposited_code/` | Verbatim copy of the Smith et al. script that defines the duration model (`2_fcwc_risk_conditional_on_differential_interrogation_length.r`), with the MIT license and README from their repository, https://github.com/thomasbryansmith/FCWC_Replication (files dated 2026-09-08). It is not run. It is included so that the inputs (their lines 284–303) and the functions (their lines 306–519) can be checked against `reply_analysis.R` line by line. |
| `figures/` | Written by the script: `Figure1_inverse_probability_problem.png`, `Figure2_fcwc_source_choice.png`, `Figure3_sensitivity_and_uncertainty.png` (300 dpi, no captions). |
| `output/results.txt` | Written by the script: everything it prints to the console. |

## Requirements

R (4.x), with the packages `dplyr` and `ragg`. Base R otherwise. Runs in about two minutes; the bootstrap is most of it.

## How to run

From the package folder:

```
Rscript reply_analysis.R
```

or open `reply_analysis.R` in RStudio, set the working directory to the source file location, and source it. Paths are relative to the package folder.

## What the script does not do

It does not re-estimate the FCWC base rate. That is a `brms` model in Smith et al.'s scripts 0 and 1 (251,000 iterations), and it reads their `data/*.xlsx` files. Nothing in the reply depends on re-running it; the reply uses their reported posterior median, 0.019, as a fixed value. Those data files are therefore not included here; they are in their repository.

## Where each number in the reply comes from

| In the reply | Value | Script section |
|---|---|---|
| Fidelity: their published posteriors at 1.53, 3.07, 6.28, 12.06, 33.96 h (0.012, 0.027, 0.10, 0.25, 0.50) | 0.0128, 0.0291, 0.1063, 0.2630, 0.5169 with π fixed | 3 |
| The False Confession Sample: fitted Redlich distribution at 6 h or more, 12 h or more | 12.1%, 3.2% | 4 |
| The False Confession Sample: Drizin & Leo Table 7, 6 h+ and 12 h+ categories | 84%, 50% (at least 50%; 12 h sits on a category boundary) | 4 |
| The Comparison Sample: weight by source type | 62.4% officer surveys, 16.0% suspect self-report, 21.6% recorded | 5 |
| The Comparison Sample: Kassin et al. (n = 601) against Kelly et al. (n = 29) | ≈ 21× | 5 |
| Thresholds Depend on the Inputs: published specification, π fixed | 10.6% at 6.28 h; crosses 10% at 6.07 h | 6 |
| Thresholds Depend on the Inputs / Figure 2: FCWC source at 12.06 h | Redlich only 26.3%; both, n-weighted 62.1%; Drizin & Leo only 72.6% | 6, 9 |
| Figure 2: the same three at 6.28 h, and their 10% crossings | 10.6%, 15.8%, 19.5%; 6.07, 5.42, 5.23 h | 6, 9 |
| Within the Same Study: Redlich true confessors as the comparison arm | 3.6% at 6.28 h; 4.0% at 12.06 h; never above ≈ 4.2%; no 10% crossing | 7 |
| Where the Data Are Not Observed: fitted non-FCWC mass above 12.06 h, above 34 h | 0.116%, 0.0024% | 8 |
| Where the Data Are Not Observed: non-FCWC dispersion +25%, +50%, −25% | crossing at 14.75 h; none (max 4.5%); 3.93 h | 8, 9 |
| Where the Data Are Not Observed: equal weights; officer surveys dropped; recorded only | 7.20 h; 8.69 h; 6.62 h | 8 |
| Where the Data Are Not Observed / Figure 3B: bootstrap 95% intervals | 5.7–14.6% at 6.28 h; 7.1–42.3% at 12.06 h | 8, 9 |

## Notes on the calculations

- **Base rate.** π = 0.019 is Smith et al.'s posterior median at full attribution (their p. 4). They average their posteriors over the base-rate distribution; the reply fixes π so that only the duration likelihoods change. That is why the fixed-π values in section 3 sit slightly above their published ones (0.1063 rather than 0.10 at 6.28 h, and a 10% crossing at 6.07 h rather than 6.28 h).
- **Crossings.** A crossing is the first duration at which the posterior rises through the threshold, found by scanning a 0.01 h grid over 0–48 h. `uniroot` is not used because the posterior is not monotone in duration for some specifications. Smith et al. interpolate on a 0.5 h grid of posterior medians.
- **Bootstrap.** Section 8 runs Smith et al.'s own `bootstrap_mixture_likelihood()` (the shading in their Fig. 2a) on both likelihoods and carries each replicate through Bayes' rule with π fixed. 400 replications; they use 1,000 for their figure. The seed and the order of calls are fixed in the script and are what produced the intervals in the text; 3.07 h is bootstrapped first, and removing it changes the random draws for the other two durations.
- **Figures.** The three PNGs in `figures/` are the figures in the submitted manuscript; the manuscript's captions are not part of the image files. The figure numbers follow the reply (Figure 2 is the FCWC-source figure, Figure 3 the sensitivity and uncertainty figure).
- **Sources for the inputs.** The study-level means, standard deviations and sample sizes are Smith et al.'s transcriptions, taken verbatim from their script. The classification of each non-FCWC row by source type (officer survey, suspect self-report, recorded interrogation) is ours, checked against each article. The Drizin & Leo (2004) counts are from their Table 7, p. 949.
