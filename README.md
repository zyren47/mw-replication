# Replication code

**Cyclical Effects of Minimum Wage Policy on Upper-Half Wage Compression in Non-Manufacturing Economies**
Ziyu Ren (University of Chicago)

This repository holds the estimation code for the paper. `paper_final_specs.R`
reproduces the main tables (1 to 4) and the appendix specifications (A1 to A3):
the headline 90/50 result, the manufacturing split, the robustness battery, the
cyclical-proxy specificity test, the channel-independence triple interaction,
the static IV appendix, the 50/10 timing diagnostic, and the full
twelve-specification wild bootstrap matrix.

## Method

Reduced-form two-way (state and year) fixed-effects regressions estimated with
`fixest::feols`, using the neighbor-state minimum wage as the policy variable
(Dube, Lester, and Reich design). Inference uses the wild cluster bootstrap
(`fwildclusterboot::boottest`, B = 9,999 Rademacher draws, null imposed),
clustered on state, given the small effective cluster counts (24 to 49).

## Files

- `paper_final_specs.R` — all regressions reported in the paper (R, using `fixest` and `fwildclusterboot`).
- `paper_final_specs.py` — a self-contained Python port producing the same specifications. The two-way fixed-effect absorption, CRV1 covariance, and restricted wild cluster bootstrap (WCR11, Rademacher, null imposed) are implemented directly in NumPy. Point estimates match the R output to numerical precision; wild bootstrap p-values match to within Monte Carlo error (about +/- 0.005 at B = 9,999), since the two languages draw independent weights.
- `config.R` — paths and bootstrap settings for the R version (`DATA_DIR`, `BOOT_B`, `BOOT_TYPE`, `BOOT_SEED`). The Python version sets the same values at the top of the script.

## Data (not included)

The scripts expect the following files in `DATA_DIR`. They are not distributed
here; all underlying sources are public.

- `panel_final.csv` — state-year panel (minimum wage, neighbor-state minimum wage, manufacturing share, and related variables).
- `hourly_wage_percentile_ratios.csv` — EPI State of Working America wage percentile ratios (50/10, 90/10, 90/50).
- `state_unemployment.csv` — annual state unemployment rates (BLS/FRED).

Sources: Economic Policy Institute (wage percentile ratios), BEA Regional
Accounts (GDP and sectoral shares), BLS/FRED (unemployment), and a hand-coded
series of state-level minimum wage overrides verified against Department of
Labor records.

## Requirements

- R: `tidyverse`, `fixest`, `fwildclusterboot`, `dqrng`.
- Python: `numpy`, `pandas` (`scipy` optional, used only for the cluster-robust t p-value).

## How to run

R: set `DATA_DIR` in `config.R`, then `source("paper_final_specs.R")`.
Python: set `DATA_DIR` at the top of `paper_final_specs.py`, then `python paper_final_specs.py`.

Each script prints every table to the console with coefficients, cluster-robust
p-values, and wild bootstrap p-values.
