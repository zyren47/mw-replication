# ==============================================================================
# paper_final_specs.R — Every Regression in the Paper (Reproducible)
# ==============================================================================
# Single source of truth for Tables 1-4 and Appendix A1-A3.
# Run AFTER 01_main_analysis_v5.R has built panel_final.csv.
# ==============================================================================

# ── Locate config.R, source it ────────────────────────────────────────────────
.config_candidates <- c("config.R", "../config.R", "./config.R")
.config_found <- FALSE
for (.c in .config_candidates) {
  if (file.exists(.c)) { source(.c); .config_found <- TRUE; break }
}
if (!.config_found) {
  stop("paper_final_specs.R: Cannot find config.R. Set working directory to ",
       "the project root (where config.R lives) before running.")
}

library(tidyverse)
library(fixest)
library(fwildclusterboot)
library(dqrng)


# ==============================================================================
# 0. BUILD ANALYSIS PANEL
# ==============================================================================

panel_raw <- read_csv(file.path(DATA_DIR, "panel_final.csv"), show_col_types = FALSE)

epi <- read_csv(file.path(DATA_DIR, "hourly_wage_percentile_ratios.csv"),
                show_col_types = FALSE)
ratios <- epi %>%
  filter(group_value %in% c("50\u201310 ratio", "50-10 ratio",
                            "90\u201310 ratio", "90-10 ratio",
                            "90\u201350 ratio", "90-50 ratio")) %>%
  mutate(metric = case_when(
    grepl("50.10", group_value) ~ "ratio_5010",
    grepl("90.10", group_value) ~ "ratio_9010",
    grepl("90.50", group_value) ~ "ratio_9050"
  )) %>%
  select(state = state_abbreviation, year, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value)

unemp <- read_csv(file.path(DATA_DIR, "state_unemployment.csv"),
                  show_col_types = FALSE)

panel <- panel_raw %>%
  left_join(ratios, by = c("state", "year")) %>%
  left_join(unemp,  by = c("state", "year")) %>%
  arrange(state, year) %>%
  group_by(state) %>%
  mutate(
    unemp_gap        = unemp_rate - mean(unemp_rate, na.rm = TRUE),
    unemp_diff       = unemp_rate - lag(unemp_rate),
    mw_x_unemp_diff  = min_wage    * unemp_diff,
    nbr_x_unemp_diff = neighbor_mw * unemp_diff,
    mw_x_unemp_gap   = min_wage    * unemp_gap,
    nbr_x_unemp_gap  = neighbor_mw * unemp_gap
  ) %>%
  ungroup()

state_manuf <- panel %>%
  group_by(state) %>%
  summarise(m_mean = mean(manuf_share, na.rm = TRUE))
high_manuf_states <- state_manuf$state[state_manuf$m_mean > median(state_manuf$m_mean)]
panel <- panel %>% mutate(high_manuf = state %in% high_manuf_states)

panel_an <- as.data.frame(droplevels(subset(panel,
                                            !state %in% c("AK","HI") &
                                              !is.na(unemp_diff) &
                                              !is.na(unemp_gap))))
panel_an$state <- as.factor(panel_an$state); panel_an$year <- as.factor(panel_an$year)

panel_low  <- panel_an[!panel_an$high_manuf, ]
panel_low$state <- droplevels(panel_low$state); panel_low$year <- droplevels(panel_low$year)

panel_high <- panel_an[panel_an$high_manuf, ]
panel_high$state <- droplevels(panel_high$state); panel_high$year <- droplevels(panel_high$year)

run_wild <- function(label, formula, data, param) {
  m  <- feols(formula, data = data, cluster = ~state)
  bt <- boottest(m, param = param, clustid = "state", B = BOOT_B,
                 type = BOOT_TYPE, impose_null = TRUE)
  cat(sprintf("%-45s N=%4d  coef=%+.5f  clu_p=%.4f  wild_p=%.4f  CI=[%+.5f, %+.5f]\n",
              label, nobs(m), coef(m)[param],
              pvalue(summary(m, cluster = ~state))[param], bt$p_val,
              bt$conf_int[1], bt$conf_int[2]))
  invisible(list(m = m, bt = bt))
}


# ==============================================================================
# TABLE 1 (MAIN) — 9050 × UnempGap, by manufacturing intensity
# ==============================================================================
cat("\n========== TABLE 1: MAIN — 9050 x UnempGap ==========\n")
set.seed(BOOT_SEED); dqset.seed(BOOT_SEED)

run_wild("Full sample",
         ratio_9050 ~ neighbor_mw + nbr_x_unemp_gap + unemp_gap | state + year,
         panel_an[!is.na(panel_an$ratio_9050), ], "nbr_x_unemp_gap")
run_wild("LOW MANUF [* main finding]",
         ratio_9050 ~ neighbor_mw + nbr_x_unemp_gap + unemp_gap | state + year,
         panel_low[!is.na(panel_low$ratio_9050), ], "nbr_x_unemp_gap")
run_wild("HIGH MANUF",
         ratio_9050 ~ neighbor_mw + nbr_x_unemp_gap + unemp_gap | state + year,
         panel_high[!is.na(panel_high$ratio_9050), ], "nbr_x_unemp_gap")


# ==============================================================================
# TABLE 2 (MAIN) — Robustness on 9050 x UnempGap x LOW MANUF
# ==============================================================================
cat("\n========== TABLE 2: ROBUSTNESS — 9050 LOW MANUF ==========\n")
set.seed(BOOT_SEED); dqset.seed(BOOT_SEED)

run_wild("Baseline (from Table 1)",
         ratio_9050 ~ neighbor_mw + nbr_x_unemp_gap + unemp_gap | state + year,
         panel_low[!is.na(panel_low$ratio_9050), ], "nbr_x_unemp_gap")

# R1 placebo
panel_low_R <- panel_low %>%
  arrange(state, year) %>% group_by(state) %>%
  mutate(unemp_gap_lead1 = lead(unemp_gap, 1),
         nbr_lead1_gap   = neighbor_mw * unemp_gap_lead1) %>%
  ungroup() %>% as.data.frame()
panel_low_R$state <- as.factor(panel_low_R$state); panel_low_R$year <- as.factor(panel_low_R$year)
run_wild("R1 Placebo (lead t+1)",
         ratio_9050 ~ neighbor_mw + nbr_lead1_gap + unemp_gap_lead1 | state + year,
         panel_low_R[!is.na(panel_low_R$ratio_9050) & !is.na(panel_low_R$nbr_lead1_gap), ],
         "nbr_lead1_gap")

# R2 exclude GFC
panel_low_nogfc <- panel_low[!panel_low$year %in% c("2008","2009","2010"), ]
panel_low_nogfc$state <- droplevels(panel_low_nogfc$state)
panel_low_nogfc$year  <- droplevels(panel_low_nogfc$year)
run_wild("R2 Exclude GFC (2008-2010)",
         ratio_9050 ~ neighbor_mw + nbr_x_unemp_gap + unemp_gap | state + year,
         panel_low_nogfc[!is.na(panel_low_nogfc$ratio_9050), ], "nbr_x_unemp_gap")

# R3 exclude CA
panel_low_noca <- panel_low[panel_low$state != "CA", ]
panel_low_noca$state <- droplevels(panel_low_noca$state)
panel_low_noca$year  <- droplevels(panel_low_noca$year)
run_wild("R3 Exclude CA",
         ratio_9050 ~ neighbor_mw + nbr_x_unemp_gap + unemp_gap | state + year,
         panel_low_noca[!is.na(panel_low_noca$ratio_9050), ], "nbr_x_unemp_gap")

# R4 state trends
run_wild("R4 + state trends",
         ratio_9050 ~ neighbor_mw + nbr_x_unemp_gap + unemp_gap + state[year] |
                       state + year,
         panel_low[!is.na(panel_low$ratio_9050), ], "nbr_x_unemp_gap")


# ==============================================================================
# TABLE 3 (MAIN) — Cyclical-proxy specificity on 9050 LOW MANUF
# ==============================================================================
cat("\n========== TABLE 3: PROXY SPECIFICITY — 9050 LOW MANUF ==========\n")
set.seed(BOOT_SEED); dqset.seed(BOOT_SEED)

run_wild("9050 x dUnemp (transitory)",
         ratio_9050 ~ neighbor_mw + nbr_x_unemp_diff + unemp_diff | state + year,
         panel_low[!is.na(panel_low$ratio_9050), ], "nbr_x_unemp_diff")
run_wild("9050 x UnempGap (persistent)",
         ratio_9050 ~ neighbor_mw + nbr_x_unemp_gap + unemp_gap | state + year,
         panel_low[!is.na(panel_low$ratio_9050), ], "nbr_x_unemp_gap")

cat("\nProxy persistence (within-state AR(1)):\n")
ar1 <- panel_an %>% group_by(state) %>% arrange(year) %>%
  summarise(gap_ar1  = cor(unemp_gap,  lag(unemp_gap),  use = "pairwise.complete.obs"),
            diff_ar1 = cor(unemp_diff, lag(unemp_diff), use = "pairwise.complete.obs"))
cat(sprintf("  UnempGap: %.3f  (persistent)\n", mean(ar1$gap_ar1,  na.rm = TRUE)))
cat(sprintf("  dUnemp:   %.3f  (transitory)\n", mean(ar1$diff_ar1, na.rm = TRUE)))


# ==============================================================================
# TABLE 4 (SUPPORTING) — Layer 3 triple interaction
# ==============================================================================
cat("\n========== TABLE 4: CHANNEL INDEPENDENCE ==========\n")

dL3 <- panel_an
dL3$manuf_dm           <- dL3$manuf_share - mean(dL3$manuf_share, na.rm = TRUE)
dL3$mw_x_manuf         <- dL3$min_wage * dL3$manuf_dm
dL3$mw_x_unemp_x_manuf <- dL3$min_wage * dL3$unemp_diff * dL3$manuf_dm

m_triple <- feols(ratio_5010 ~ min_wage + unemp_diff + manuf_dm +
                                mw_x_unemp_diff + mw_x_manuf + mw_x_unemp_x_manuf |
                                state + year,
                  data = dL3[!is.na(dL3$ratio_5010), ], cluster = ~state)
print(summary(m_triple, cluster = ~state))
# Key: mw_x_unemp_x_manuf p ≈ 0.67 (null) → channels independent


# ==============================================================================
# APPENDIX A1 — Layer 1 static IV (wild fragile)
# ==============================================================================
cat("\n========== APPENDIX A1: STATIC IV REPLICATION ==========\n")
set.seed(BOOT_SEED); dqset.seed(BOOT_SEED)

dA1 <- as.data.frame(droplevels(subset(panel, !state %in% c("AK","HI"))))
dA1$state <- as.factor(dA1$state); dA1$year <- as.factor(dA1$year)

cat(sprintf("\n%-12s %5s %8s %8s %8s %8s %8s\n",
            "outcome","N","IV_coef","IV_clu_p","RF_coef","RF_clu_p","RF_wild_p"))
cat(strrep("-", 70), "\n")
for (out in c("ratio_5010","ratio_9050","ratio_9010")) {
  d_sub <- dA1[!is.na(dA1[[out]]), ]
  f_iv  <- as.formula(sprintf("%s ~ 1 | state + year | min_wage ~ neighbor_mw", out))
  m_iv  <- feols(f_iv, data = d_sub, cluster = ~state)
  iv_c  <- coef(m_iv)["fit_min_wage"]
  iv_p  <- pvalue(summary(m_iv, cluster = ~state))["fit_min_wage"]
  f_rf  <- as.formula(sprintf("%s ~ neighbor_mw | state + year", out))
  m_rf  <- feols(f_rf, data = d_sub, cluster = ~state)
  rf_c  <- coef(m_rf)["neighbor_mw"]
  rf_p  <- pvalue(summary(m_rf, cluster = ~state))["neighbor_mw"]
  bt    <- boottest(m_rf, param = "neighbor_mw", clustid = "state",
                    B = BOOT_B, type = BOOT_TYPE, impose_null = TRUE)
  cat(sprintf("%-12s %5d %+8.4f %8.4f %+8.4f %8.4f %8.4f\n",
              out, nobs(m_rf), iv_c, iv_p, rf_c, rf_p, bt$p_val))
}


# ==============================================================================
# APPENDIX A2 — 5010 cyclical with placebo failure disclosure
# ==============================================================================
cat("\n========== APPENDIX A2: 5010 TIMING DIAGNOSTIC ==========\n")

dA2 <- panel_low %>%
  arrange(state, year) %>% group_by(state) %>%
  mutate(unemp_diff_lead1 = lead(unemp_diff, 1),
         unemp_diff_lead2 = lead(unemp_diff, 2),
         unemp_diff_lag1  = lag(unemp_diff,  1),
         nbr_lead1 = neighbor_mw * unemp_diff_lead1,
         nbr_lead2 = neighbor_mw * unemp_diff_lead2,
         nbr_lag1  = neighbor_mw * unemp_diff_lag1) %>%
  ungroup() %>% as.data.frame()
dA2$state <- as.factor(dA2$state); dA2$year <- as.factor(dA2$year)

set.seed(BOOT_SEED); dqset.seed(BOOT_SEED)

run_wild("Contemp (t)",
         ratio_5010 ~ neighbor_mw + nbr_x_unemp_diff + unemp_diff | state + year,
         dA2[!is.na(dA2$ratio_5010), ], "nbr_x_unemp_diff")
run_wild("Lead (t+1) - placebo",
         ratio_5010 ~ neighbor_mw + nbr_lead1 + unemp_diff_lead1 | state + year,
         dA2[!is.na(dA2$ratio_5010) & !is.na(dA2$nbr_lead1), ], "nbr_lead1")
run_wild("Lead (t+2) - far placebo",
         ratio_5010 ~ neighbor_mw + nbr_lead2 + unemp_diff_lead2 | state + year,
         dA2[!is.na(dA2$ratio_5010) & !is.na(dA2$nbr_lead2), ], "nbr_lead2")
run_wild("Lag (t-1)",
         ratio_5010 ~ neighbor_mw + nbr_lag1 + unemp_diff_lag1 | state + year,
         dA2[!is.na(dA2$ratio_5010) & !is.na(dA2$nbr_lag1), ], "nbr_lag1")

cat("\n--- Horse race: contemp + lead(t+1) jointly ---\n")
m_horse <- feols(ratio_5010 ~ neighbor_mw + nbr_x_unemp_diff + unemp_diff +
                                 nbr_lead1 + unemp_diff_lead1 | state + year,
                 data = dA2[!is.na(dA2$ratio_5010) & !is.na(dA2$nbr_lead1), ],
                 cluster = ~state)
print(summary(m_horse, cluster = ~state))
bt_contemp <- boottest(m_horse, param = "nbr_x_unemp_diff", clustid = "state",
                       B = BOOT_B, type = BOOT_TYPE, impose_null = TRUE)
bt_lead    <- boottest(m_horse, param = "nbr_lead1", clustid = "state",
                       B = BOOT_B, type = BOOT_TYPE, impose_null = TRUE)
cat(sprintf("\nHorse race wild p: contemp=%.4f, lead=%.4f\n",
            bt_contemp$p_val, bt_lead$p_val))


# ==============================================================================
# APPENDIX A3 — Full 12-spec Layer 2 wild bootstrap matrix
# ==============================================================================
cat("\n========== APPENDIX A3: FULL 12-SPEC MATRIX ==========\n")
set.seed(BOOT_SEED); dqset.seed(BOOT_SEED)

specs <- expand.grid(
  outcome = c("ratio_5010","ratio_9050"),
  proxy   = c("dUnemp","UnempGap"),
  sample  = c("Full","Low","High"),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(specs))) {
  o  <- specs$outcome[i]; p <- specs$proxy[i]; s <- specs$sample[i]
  d_use <- switch(s, "Full" = panel_an, "Low" = panel_low, "High" = panel_high)
  d_use <- d_use[!is.na(d_use[[o]]), ]
  if (p == "dUnemp") {
    int_var <- "nbr_x_unemp_diff"; cyc_main <- "unemp_diff"
  } else {
    int_var <- "nbr_x_unemp_gap";  cyc_main <- "unemp_gap"
  }
  f <- as.formula(sprintf("%s ~ neighbor_mw + %s + %s | state + year",
                          o, int_var, cyc_main))
  run_wild(sprintf("%s x %s x %s", o, p, s), f, d_use, int_var)
}


cat("\n========== PAPER SPECIFICATIONS COMPLETE ==========\n\n")
cat("Expected numbers (verified 2026-05-23):\n")
cat("  Table 1 (Main):     9050 x UnempGap x {Full=0.046, Low=0.041, High=0.665}\n")
cat("  Table 2 (Robust):   R1=0.293, R2=0.033, R3=0.058, R4=0.038\n")
cat("  Table 3 (Proxy):    dUnemp=0.797 (ns), UnempGap=0.041 (sig)\n")
cat("  Table 4 (Triple):   mw_x_unemp_x_manuf p=0.67\n")
cat("  App A1 (L1):        all wild p > 0.10 (wild fragile)\n")
cat("  App A2 (5010):      contemp 0.006, lead-t+1 0.012, lead-t+2 0.033\n")
cat("  App A3 (12-grid):   4 robust cells, all in {Low manuf, Full}\n")
