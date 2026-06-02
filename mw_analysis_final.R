# ==============================================================================
# Cyclical Effects of Minimum Wage Policy on Employment and Inequality
# Ziyu Ren | University of Chicago, MA Economics
#
# Replication code: panel construction, TWFE estimation, robustness checks,
# neighbor-state IV, heterogeneity analysis, and tradeoff frontier plot.
#
# Data (place in DATA_DIR):
#   SAGDP1__ALL_AREAS_1997_2024.csv   BEA real GDP by state
#   SAGDP2__ALL_AREAS_1997_2024.csv   BEA GDP by industry
#   Frank_StatesUS_2022.xls           Top 10% income share (Frank 2022)
#   hourly_wage_percentile_ratios.csv EPI 50-10 wage ratio
#
# Packages: tidyverse, fixest, readxl
# ==============================================================================

library(tidyverse)
library(fixest)
library(readxl)

DATA_DIR <- "/Users/ProjectPei/Desktop/"

# ==============================================================================
# 0. CONSTANTS
# ==============================================================================

FEDERAL_MW <- tibble(
  year   = 1997:2018,
  fed_mw = c(4.75, 5.15, 5.15, 5.15, 5.15, 5.15, 5.15, 5.15, 5.15, 5.15,
             5.85, 6.55, 7.25, 7.25, 7.25, 7.25, 7.25, 7.25, 7.25, 7.25,
             7.25, 7.25)
)

# 224 state-year cells where effective state MW exceeds the federal floor.
# Cross-checked against DOL Wage and Hour Division records and
# state ballot-measure implementation dates.
STATE_MW_OVERRIDES <- tribble(
  ~state, ~year, ~state_mw,
  "CA",2001,6.25, "CA",2002,6.75, "CA",2003,6.75, "CA",2004,6.75,
  "CA",2005,6.75, "CA",2006,6.75, "CA",2007,7.50, "CA",2008,8.00,
  "CA",2009,8.00, "CA",2010,8.00, "CA",2011,8.00, "CA",2012,8.00,
  "CA",2013,8.00, "CA",2014,9.00, "CA",2015,9.00, "CA",2016,10.00,
  "CA",2017,10.50,"CA",2018,11.00,
  "NY",2005,6.00, "NY",2006,6.75, "NY",2007,7.15, "NY",2008,7.15,
  "NY",2014,8.00, "NY",2015,8.75, "NY",2016,9.00, "NY",2017,9.70,
  "NY",2018,10.40,
  "WA",1999,5.70, "WA",2000,6.50, "WA",2001,6.72, "WA",2002,6.90,
  "WA",2003,7.01, "WA",2004,7.16, "WA",2005,7.35, "WA",2006,7.63,
  "WA",2007,7.93, "WA",2008,8.07, "WA",2009,8.55, "WA",2010,8.55,
  "WA",2011,8.67, "WA",2012,9.04, "WA",2013,9.19, "WA",2014,9.32,
  "WA",2015,9.47, "WA",2016,9.47, "WA",2017,11.00,"WA",2018,11.50,
  "MA",2007,7.50, "MA",2008,8.00, "MA",2015,9.00, "MA",2016,10.00,
  "MA",2017,11.00,"MA",2018,11.00,
  "OR",1999,6.00, "OR",2000,6.50, "OR",2001,6.50, "OR",2002,6.50,
  "OR",2003,6.90, "OR",2004,7.05, "OR",2005,7.25, "OR",2006,7.50,
  "OR",2007,7.80, "OR",2008,7.95, "OR",2009,8.40, "OR",2010,8.40,
  "OR",2011,8.50, "OR",2012,8.80, "OR",2013,8.95, "OR",2014,9.10,
  "OR",2015,9.25, "OR",2016,9.25, "OR",2017,9.75, "OR",2018,10.25,
  "IL",2004,5.50, "IL",2005,6.50, "IL",2006,6.50, "IL",2007,7.50,
  "IL",2008,7.75, "IL",2009,8.00, "IL",2010,8.25, "IL",2011,8.25,
  "IL",2012,8.25, "IL",2013,8.25, "IL",2014,8.25, "IL",2015,8.25,
  "IL",2016,8.25, "IL",2017,8.25, "IL",2018,8.25,
  "CT",2014,8.70, "CT",2015,9.15, "CT",2016,9.60, "CT",2017,10.10,
  "CT",2018,10.10,
  "MN",2014,8.00, "MN",2015,9.00, "MN",2016,9.00, "MN",2017,9.50,
  "MN",2018,9.65,
  "NJ",2014,8.25, "NJ",2015,8.38, "NJ",2016,8.38, "NJ",2017,8.44,
  "NJ",2018,8.60,
  "CO",2007,6.85, "CO",2008,7.02, "CO",2009,7.28, "CO",2010,7.24,
  "CO",2011,7.36, "CO",2012,7.64, "CO",2013,7.78, "CO",2014,8.00,
  "CO",2015,8.23, "CO",2016,8.31, "CO",2017,9.30, "CO",2018,10.20,
  "AZ",2007,6.75, "AZ",2008,6.90, "AZ",2009,7.25, "AZ",2010,7.25,
  "AZ",2011,7.35, "AZ",2012,7.65, "AZ",2013,7.80, "AZ",2014,7.90,
  "AZ",2015,8.05, "AZ",2016,8.05, "AZ",2017,10.00,"AZ",2018,10.50,
  "VT",2007,7.53, "VT",2008,7.68, "VT",2009,8.06, "VT",2010,8.06,
  "VT",2011,8.15, "VT",2012,8.46, "VT",2013,8.60, "VT",2014,8.73,
  "VT",2015,9.15, "VT",2016,9.60, "VT",2017,10.00,"VT",2018,10.50,
  "MI",2006,6.95, "MI",2007,7.15, "MI",2008,7.40, "MI",2014,8.15,
  "MI",2015,8.15, "MI",2016,8.50, "MI",2017,8.90, "MI",2018,9.25,
  "MD",2015,8.00, "MD",2016,8.75, "MD",2017,9.25, "MD",2018,9.25,
  "NE",2015,8.00, "NE",2016,9.00, "NE",2017,9.00, "NE",2018,9.00,
  "SD",2015,8.50, "SD",2016,8.55, "SD",2017,8.65, "SD",2018,8.85,
  "AK",2003,7.15, "AK",2014,8.75, "AK",2015,9.75, "AK",2016,9.75,
  "AK",2017,9.80, "AK",2018,9.84,
  "DC",1999,6.15, "DC",2003,6.60, "DC",2006,7.00, "DC",2008,7.55,
  "DC",2009,8.25, "DC",2010,8.25, "DC",2011,8.25, "DC",2012,8.25,
  "DC",2013,8.25, "DC",2014,9.50, "DC",2015,10.50,"DC",2016,11.50,
  "DC",2017,12.50,"DC",2018,13.25,
  "RI",2014,8.00, "RI",2015,9.00, "RI",2016,9.60, "RI",2017,9.60,
  "RI",2018,10.10,
  "HI",2003,6.25, "HI",2004,6.25, "HI",2005,6.75, "HI",2006,7.25,
  "HI",2007,7.25, "HI",2015,7.75, "HI",2016,8.50, "HI",2017,9.25,
  "HI",2018,10.10,
  "NM",2007,5.15, "NM",2008,6.50, "NM",2009,7.50,
  "FL",2005,6.15, "FL",2006,6.40, "FL",2007,6.67, "FL",2008,6.79,
  "FL",2009,7.21, "FL",2011,7.31, "FL",2012,7.67, "FL",2013,7.79,
  "FL",2014,7.93, "FL",2015,8.05, "FL",2016,8.05, "FL",2017,8.10,
  "FL",2018,8.25,
  "MT",2007,6.15, "MT",2008,6.25, "MT",2009,6.90, "MT",2015,8.05,
  "MT",2016,8.05, "MT",2017,8.15, "MT",2018,8.30,
  "NV",2009,7.55, "NV",2010,8.25, "NV",2011,8.25, "NV",2015,8.25,
  "MO",2007,6.65, "MO",2008,6.65, "MO",2009,7.05,
  "OH",2006,6.85, "OH",2007,7.00, "OH",2008,7.00, "OH",2009,7.30,
  "OH",2010,7.30, "OH",2011,7.40, "OH",2012,7.70, "OH",2013,7.85,
  "OH",2014,8.00, "OH",2015,8.10, "OH",2016,8.15, "OH",2017,8.15,
  "OH",2018,8.30,
  "ME",2017,9.00, "ME",2018,10.00,
  "WV",2016,8.75,
  "KS",2010,7.25
)

# Neighbor-state adjacency for IV construction
ADJACENCY <- list(
  AL=c("FL","GA","MS","TN"),    AK=character(0),
  AZ=c("CA","CO","NM","NV","UT"),
  AR=c("LA","MO","MS","OK","TN","TX"), CA=c("AZ","NV","OR"),
  CO=c("AZ","KS","NE","NM","OK","UT","WY"), CT=c("MA","NY","RI"),
  DE=c("MD","NJ","PA"),         FL=c("AL","GA"),
  GA=c("AL","FL","NC","SC","TN"),         HI=character(0),
  ID=c("MT","NV","OR","UT","WA","WY"),    IL=c("IN","IA","KY","MO","WI"),
  IN=c("IL","KY","MI","OH"),    IA=c("IL","MN","MO","NE","SD","WI"),
  KS=c("CO","MO","NE","OK"),    KY=c("IL","IN","MO","OH","TN","VA","WV"),
  LA=c("AR","MS","TX"),         ME=c("NH"),
  MD=c("DE","PA","VA","WV","DC"),         MA=c("CT","NH","NY","RI","VT"),
  MI=c("IN","OH","WI"),         MN=c("IA","ND","SD","WI"),
  MS=c("AL","AR","LA","TN"),    MO=c("AR","IL","IA","KS","KY","NE","OK","TN"),
  MT=c("ID","ND","SD","WY"),    NE=c("CO","IA","KS","MO","SD","WY"),
  NV=c("AZ","CA","ID","OR","UT"),         NH=c("MA","ME","VT"),
  NJ=c("DE","NY","PA"),         NM=c("AZ","CO","OK","TX"),
  NY=c("CT","MA","NJ","PA","VT"),         NC=c("GA","SC","TN","VA"),
  ND=c("MN","MT","SD"),         OH=c("IN","KY","MI","PA","WV"),
  OK=c("AR","CO","KS","MO","NM","TX"),    OR=c("CA","ID","NV","WA"),
  PA=c("DE","MD","NJ","NY","OH","WV"),    RI=c("CT","MA"),
  SC=c("GA","NC"),              SD=c("IA","MN","MT","ND","NE","WY"),
  TN=c("AL","AR","GA","KY","MS","MO","NC","VA"), TX=c("AR","LA","NM","OK"),
  UT=c("AZ","CO","ID","NV","NM","WY"),    VT=c("MA","NH","NY"),
  VA=c("KY","MD","NC","TN","WV","DC"),    WA=c("ID","OR"),
  WV=c("KY","MD","OH","PA","VA"),         WI=c("IL","IA","MI","MN"),
  WY=c("CO","ID","MT","NE","SD","UT"),    DC=c("MD","VA")
)

FIPS_TO_STATE <- c(
  "01000"="AL","02000"="AK","04000"="AZ","05000"="AR","06000"="CA",
  "08000"="CO","09000"="CT","10000"="DE","12000"="FL","13000"="GA",
  "15000"="HI","16000"="ID","17000"="IL","18000"="IN","19000"="IA",
  "20000"="KS","21000"="KY","22000"="LA","23000"="ME","24000"="MD",
  "25000"="MA","26000"="MI","27000"="MN","28000"="MS","29000"="MO",
  "30000"="MT","31000"="NE","32000"="NV","33000"="NH","34000"="NJ",
  "35000"="NM","36000"="NY","37000"="NC","38000"="ND","39000"="OH",
  "40000"="OK","41000"="OR","42000"="PA","44000"="RI","45000"="SC",
  "46000"="SD","47000"="TN","48000"="TX","49000"="UT","50000"="VT",
  "51000"="VA","53000"="WA","54000"="WV","55000"="WI","56000"="WY",
  "11000"="DC"
)

NAME_TO_ABBR <- c(
  "Alabama"="AL","Alaska"="AK","Arizona"="AZ","Arkansas"="AR",
  "California"="CA","Colorado"="CO","Connecticut"="CT","Delaware"="DE",
  "Florida"="FL","Georgia"="GA","Hawaii"="HI","Idaho"="ID",
  "Illinois"="IL","Indiana"="IN","Iowa"="IA","Kansas"="KS",
  "Kentucky"="KY","Louisiana"="LA","Maine"="ME","Maryland"="MD",
  "Massachusetts"="MA","Michigan"="MI","Minnesota"="MN","Mississippi"="MS",
  "Missouri"="MO","Montana"="MT","Nebraska"="NE","Nevada"="NV",
  "New Hampshire"="NH","New Jersey"="NJ","New Mexico"="NM","New York"="NY",
  "North Carolina"="NC","North Dakota"="ND","Ohio"="OH","Oklahoma"="OK",
  "Oregon"="OR","Pennsylvania"="PA","Rhode Island"="RI",
  "South Carolina"="SC","South Dakota"="SD","Tennessee"="TN","Texas"="TX",
  "Utah"="UT","Vermont"="VT","Virginia"="VA","Washington"="WA",
  "West Virginia"="WV","Wisconsin"="WI","Wyoming"="WY",
  "District of Columbia"="DC"
)

YR_COLS <- as.character(1997:2018)

# ==============================================================================
# 1. DATA LOADING & PANEL CONSTRUCTION
# ==============================================================================

load_gdp <- function() {
  bea1 <- read_csv(paste0(DATA_DIR, "SAGDP1__ALL_AREAS_1997_2024.csv"),
                   show_col_types = FALSE)
  gdp <- bea1 %>%
    filter(GeoFIPS %in% names(FIPS_TO_STATE), LineCode == 1) %>%
    mutate(state = FIPS_TO_STATE[GeoFIPS]) %>%
    select(state, all_of(YR_COLS)) %>%
    pivot_longer(-state, names_to = "year", values_to = "real_gdp") %>%
    mutate(year     = as.integer(year),
           real_gdp = suppressWarnings(as.numeric(real_gdp))) %>%
    filter(year >= 1998) %>%
    arrange(state, year) %>%
    group_by(state) %>%
    mutate(gdp_growth = (real_gdp - lag(real_gdp)) / lag(real_gdp)) %>%
    ungroup()
  gdp
}

load_sector_shares <- function() {
  bea2 <- read_csv(paste0(DATA_DIR, "SAGDP2__ALL_AREAS_1997_2024.csv"),
                   show_col_types = FALSE) %>%
    mutate(Description = str_trim(Description)) %>%
    filter(GeoFIPS %in% names(FIPS_TO_STATE)) %>%
    mutate(state = FIPS_TO_STATE[GeoFIPS])

  pivot_sector <- function(desc_vec, colname) {
    bea2 %>%
      filter(Description %in% desc_vec) %>%
      select(state, all_of(YR_COLS)) %>%
      group_by(state) %>%
      summarise(across(everything(),
                       ~ suppressWarnings(sum(as.numeric(.), na.rm = TRUE))),
                .groups = "drop") %>%
      pivot_longer(-state, names_to = "year", values_to = colname) %>%
      mutate(year = as.integer(year))
  }

  total <- pivot_sector("All industry total",    "gdp_total")
  fin   <- pivot_sector("Finance and insurance", "gdp_fin")
  manuf <- pivot_sector(c("Durable goods manufacturing",
                           "Nondurable goods manufacturing"), "gdp_manuf")

  total %>%
    left_join(fin,   by = c("state","year")) %>%
    left_join(manuf, by = c("state","year")) %>%
    mutate(finance_share = gdp_fin   / gdp_total * 100,
           manuf_share   = gdp_manuf / gdp_total * 100) %>%
    filter(year >= 1998) %>%
    select(state, year, finance_share, manuf_share)
}

load_frank <- function() {
  frank <- read_excel(paste0(DATA_DIR, "Frank_StatesUS_2022.xls"),
                      sheet = "Top Income ShareUS", col_names = FALSE)
  colnames(frank) <- c("year","st","state_name","top10us","top5us","top1us",
                       "top05us","top01us","x","xx")
  frank %>%
    mutate(state   = NAME_TO_ABBR[state_name],
           year    = as.integer(year),
           top10us = suppressWarnings(as.numeric(top10us))) %>%
    filter(year >= 1998, year <= 2018, !is.na(state)) %>%
    select(state, year, top10us)
}

load_wage_ratio <- function() {
  epi <- read_csv(paste0(DATA_DIR, "hourly_wage_percentile_ratios.csv"),
                  show_col_types = FALSE)
  epi %>%
    filter(geo_type == "state",
           str_detect(group_value, regex("50.10 ratio", ignore_case = TRUE)),
           year >= 1998, year <= 2018) %>%
    mutate(state      = NAME_TO_ABBR[geo_name],
           ratio_5010 = suppressWarnings(as.numeric(value))) %>%
    filter(!is.na(state)) %>%
    select(state, year, ratio_5010)
}

build_minimum_wage <- function() {
  expand_grid(state = names(ADJACENCY), year = 1998:2018) %>%
    left_join(FEDERAL_MW,         by = "year") %>%
    left_join(STATE_MW_OVERRIDES, by = c("state","year")) %>%
    mutate(min_wage = pmax(fed_mw, replace_na(state_mw, 0))) %>%
    select(state, year, min_wage)
}

build_neighbor_mw <- function(mw_df) {
  # Instrument: average MW across contiguous states.
  # Exclusion restriction: neighbor legislatures set wages through independent
  # political processes uncorrelated with own-state income distribution.
  mw_lookup <- mw_df %>% select(state, year, min_wage)
  mw_df %>%
    rowwise() %>%
    mutate(neighbor_mw = {
      nbrs <- ADJACENCY[[state]]
      if (length(nbrs) == 0) {
        NA_real_
      } else {
        vals <- mw_lookup$min_wage[mw_lookup$state %in% nbrs &
                                     mw_lookup$year == year]
        if (length(vals) == 0) NA_real_ else mean(vals, na.rm = TRUE)
      }
    }) %>%
    ungroup() %>%
    select(state, year, neighbor_mw)
}

build_panel <- function() {
  gdp     <- load_gdp()
  sectors <- load_sector_shares()
  frank   <- load_frank()
  wr      <- load_wage_ratio()
  mw      <- build_minimum_wage()
  nbr     <- build_neighbor_mw(mw)

  panel <- gdp %>%
    left_join(mw,      by = c("state","year")) %>%
    left_join(sectors, by = c("state","year")) %>%
    left_join(frank,   by = c("state","year")) %>%
    left_join(wr,      by = c("state","year")) %>%
    left_join(nbr,     by = c("state","year")) %>%
    filter(!is.na(gdp_growth), !is.na(min_wage)) %>%
    arrange(state, year) %>%
    group_by(state) %>%
    mutate(mw_lag1   = lag(min_wage),
           mwlag_gdp = mw_lag1 * gdp_growth) %>%
    ungroup() %>%
    mutate(
      mw_gdp          = min_wage * gdp_growth,
      gdp_growth_sq   = gdp_growth ^ 2,
      neighbor_mw_gdp = neighbor_mw * gdp_growth,
      finance_gdp     = finance_share * gdp_growth,
      recession       = as.integer(year %in% c(2008, 2009)),
      mw_rec          = min_wage * recession
    ) %>%
    group_by(year) %>%
    mutate(mw_yr_med = median(min_wage, na.rm = TRUE)) %>%
    ungroup() %>%
    # Kaitz index: captures whether MW is actually binding in local labor markets
    mutate(mw_relative = min_wage / mw_yr_med)

  cat(sprintf("Panel: %d obs | %d states | %d-%d\n",
              nrow(panel), n_distinct(panel$state),
              min(panel$year, na.rm = TRUE), max(panel$year, na.rm = TRUE)))
  cat(sprintf("  top10us:    %d obs\n", sum(!is.na(panel$top10us))))
  cat(sprintf("  ratio_5010: %d obs\n", sum(!is.na(panel$ratio_5010))))
  panel
}

# ==============================================================================
# 2. SUMMARY STATISTICS
# ==============================================================================

summary_stats <- function(panel) {
  cat("\n", strrep("=", 68), "\n")
  cat("Summary Statistics\n")
  cat(strrep("=", 68), "\n\n")

  vars <- c(
    min_wage      = "Effective minimum wage ($)",
    gdp_growth    = "Real GDP growth rate",
    top10us       = "Top 10% income share",
    ratio_5010    = "50-10 wage ratio",
    finance_share = "Finance share of GDP (%)",
    manuf_share   = "Manufacturing share of GDP (%)",
    mw_relative   = "Kaitz index"
  )

  map_dfr(names(vars), function(v) {
    x <- panel[[v]]
    tibble(
      Variable = vars[[v]],
      N        = sum(!is.na(x)),
      Mean     = round(mean(x, na.rm = TRUE), 3),
      SD       = round(sd(x,   na.rm = TRUE), 3),
      Min      = round(min(x,  na.rm = TRUE), 3),
      Median   = round(median(x, na.rm = TRUE), 3),
      Max      = round(max(x,  na.rm = TRUE), 3)
    )
  }) %>% print(n = 10)
}

# ==============================================================================
# 3. ANALYSIS
# ==============================================================================

run_layer1 <- function(panel) {
  cat("\n", strrep("=", 68), "\n")
  cat("Layer 1: Cyclical Asymmetry — Top 10% Income Share\n")
  cat(strrep("=", 68), "\n\n")

  p <- panel %>% filter(!is.na(top10us))

  m1 <- feols(top10us ~ min_wage + gdp_growth          | state + year,
              data = p, cluster = ~state)
  m2 <- feols(top10us ~ min_wage + gdp_growth + mw_gdp | state + year,
              data = p, cluster = ~state)

  print(etable(m1, m2,
               keep   = c("min_wage","gdp_growth","mw_gdp"),
               digits = 4,
               title  = "Top 10% Income Share"))

  b0 <- coef(m2)["min_wage"]
  b1 <- coef(m2)["mw_gdp"]
  cat(sprintf("\nMW×GDP = %.4f (SE=%.4f) | threshold: GDP growth < %.1f%%\n",
              b1, se(m2)["mw_gdp"], (-b0/b1) * 100))

  invisible(list(basic = m1, interact = m2))
}

run_robustness <- function(panel) {
  cat("\n", strrep("=", 68), "\n")
  cat("Robustness: MW×GDP Across Six Specifications\n")
  cat(strrep("=", 68), "\n\n")

  p     <- panel %>% filter(!is.na(top10us))
  p_lag <- p     %>% filter(!is.na(mw_lag1))

  models <- list(
    "Baseline"       = feols(top10us ~ min_wage + gdp_growth + mw_gdp
                              | state + year, data = p, cluster = ~state),
    "+ GDP²"         = feols(top10us ~ min_wage + gdp_growth + mw_gdp +
                               gdp_growth_sq
                              | state + year, data = p, cluster = ~state),
    "+ Finance"      = feols(top10us ~ min_wage + gdp_growth + mw_gdp +
                               finance_share
                              | state + year, data = p, cluster = ~state),
    "+ Finance×GDP"  = feols(top10us ~ min_wage + gdp_growth + mw_gdp +
                               finance_share + finance_gdp
                              | state + year, data = p, cluster = ~state),
    "+ GDP²+Finance" = feols(top10us ~ min_wage + gdp_growth + mw_gdp +
                               gdp_growth_sq + finance_share
                              | state + year, data = p, cluster = ~state),
    "Lagged MW"      = feols(top10us ~ mw_lag1 + gdp_growth + mwlag_gdp +
                               gdp_growth_sq + finance_share
                              | state + year, data = p_lag, cluster = ~state)
  )

  keys <- c(rep("mw_gdp", 5), "mwlag_gdp")
  results <- map2_dfr(models, keys, function(m, k) {
    tibble(coef = coef(m)[k], se = se(m)[k], pval = pvalue(m)[k])
  }) %>%
    mutate(spec = names(models),
           sig  = case_when(pval < .01 ~ "***",
                            pval < .05 ~ "**",
                            pval < .10 ~ "*",
                            TRUE       ~ "(ns)")) %>%
    select(spec, coef, se, pval, sig)

  print(results)
  cat(sprintf("\nRange: %.4f – %.4f | All p < 0.01: %s\n",
              min(results$coef, na.rm = TRUE),
              max(results$coef, na.rm = TRUE),
              if (all(results$pval < .01)) "YES" else "NO"))
  invisible(models)
}

run_iv <- function(panel, outcome = "top10us") {
  cat("\n", strrep("=", 68), "\n")
  cat(sprintf("IV: Neighbor-State MW — outcome: %s\n", outcome))
  cat(strrep("=", 68), "\n\n")

  p <- panel %>% filter(!is.na(.data[[outcome]]), !is.na(neighbor_mw))
  cat(sprintf("Sample: %d obs, %d states (AK and HI excluded)\n\n",
              nrow(p), n_distinct(p$state)))

  fml_ols <- as.formula(sprintf(
    "%s ~ min_wage + mw_gdp + gdp_growth | state + year", outcome))
  fml_iv  <- as.formula(sprintf(
    "%s ~ gdp_growth | state + year | min_wage + mw_gdp ~ neighbor_mw + neighbor_mw_gdp",
    outcome))

  m_ols <- feols(fml_ols, data = p, cluster = ~state)
  m_iv  <- feols(fml_iv,  data = p, cluster = ~state)

  cat("First-stage F-statistics:\n")
  print(fitstat(m_iv, "ivf"))
  cat("\nOLS vs IV:\n")
  print(etable(m_ols, m_iv,
               keep   = c("min_wage","mw_gdp","fit_min_wage","fit_mw_gdp"),
               digits = 4))
  invisible(list(ols = m_ols, iv = m_iv))
}

run_heterogeneity <- function(panel, outcome = "top10us") {
  cat("\n", strrep("=", 68), "\n")
  cat(sprintf("Heterogeneity Analysis — outcome: %s\n", outcome))
  cat(strrep("=", 68), "\n")

  p <- panel %>% filter(!is.na(.data[[outcome]]))

  state_chars <- p %>%
    group_by(state) %>%
    summarise(kaitz_mean = mean(mw_relative, na.rm = TRUE),
              gdp_vol    = sd(gdp_growth,    na.rm = TRUE),
              manuf_mean = mean(manuf_share,  na.rm = TRUE),
              .groups    = "drop") %>%
    mutate(high_kaitz = kaitz_mean > median(kaitz_mean, na.rm = TRUE),
           high_vol   = gdp_vol    > median(gdp_vol,    na.rm = TRUE),
           high_manuf = manuf_mean > median(manuf_mean, na.rm = TRUE))

  p <- p %>% left_join(state_chars, by = "state")
  fml <- as.formula(sprintf(
    "%s ~ min_wage + gdp_growth + mw_gdp + gdp_growth_sq | state + year", outcome))

  dims <- list(
    list(var="high_kaitz", label="Kaitz index",
         pred="High Kaitz → asymmetry stronger where MW is binding"),
    list(var="high_vol",   label="GDP volatility",
         pred="More volatile → larger cycle shocks"),
    list(var="high_manuf", label="Manufacturing share",
         pred="Higher labor demand elasticity → cost more cycle-sensitive")
  )

  for (d in dims) {
    cat(sprintf("\n  %s (%s)\n", d$label, d$pred))
    for (grp in c(TRUE, FALSE)) {
      sub <- p %>% filter(.data[[d$var]] == grp)
      m   <- feols(fml, data = sub, cluster = ~state)
      b   <- coef(m)["mw_gdp"]
      pv  <- pvalue(m)["mw_gdp"]
      sig <- if (pv<.01)"***" else if (pv<.05)"**" else if (pv<.10)"*" else "(ns)"
      cat(sprintf("    %-4s (N=%5d): MW×GDP = %7.4f%-4s p=%.3f\n",
                  if (grp) "High" else "Low", nrow(sub), b, sig, pv))
    }
  }
}

plot_tradeoff <- function(panel) {
  p  <- panel %>% filter(!is.na(top10us))
  m  <- feols(top10us ~ min_wage + gdp_growth + mw_gdp | state + year,
              data = p, cluster = ~state)
  b0 <- coef(m)["min_wage"]
  b1 <- coef(m)["mw_gdp"]

  tibble(g = seq(-0.06, 0.08, length.out = 300), me = b0 + b1 * g) %>%
    ggplot(aes(g * 100, me)) +
    annotate("rect", xmin=-6, xmax=0, ymin=-Inf, ymax=Inf, fill="red",   alpha=0.07) +
    annotate("rect", xmin= 0, xmax=8, ymin=-Inf, ymax=Inf, fill="green", alpha=0.07) +
    geom_hline(yintercept=0, linetype="dashed", color="gray50", linewidth=0.7) +
    geom_vline(xintercept=(-b0/b1)*100, linetype="dotted",
               color="#c0392b", linewidth=0.9) +
    geom_line(color="#c0392b", linewidth=2) +
    annotate("text", x=-4.5, y=max(b0+b1*seq(-.06,.08,.001))*0.85,
             label="Contraction", color="#c0392b", fontstyle="italic", size=4) +
    annotate("text", x= 5.5, y=max(b0+b1*seq(-.06,.08,.001))*0.85,
             label="Expansion",   color="darkgreen", fontstyle="italic", size=4) +
    labs(x       = "GDP Growth Rate (%)",
         y       = expression(partialdiff*"Top10%Share" / partialdiff*"MW"),
         title   = "Equity-Efficiency Tradeoff Frontier",
         subtitle= sprintf("Threshold: GDP growth < %.1f%%", (-b0/b1)*100),
         caption = "50-state panel, 1998-2018. TWFE, SE clustered by state.") +
    theme_bw(base_size=13) +
    theme(plot.title   = element_text(face="bold", size=14),
          plot.subtitle= element_text(size=9, color="gray40"),
          plot.caption = element_text(size=8, color="gray50"))
}

# ==============================================================================
# MAIN
# ==============================================================================

panel <- build_panel()
write_csv(panel, paste0(DATA_DIR, "panel_final.csv"))

summary_stats(panel)
r1  <- run_layer1(panel)
rob <- run_robustness(panel)
iv1 <- run_iv(panel, "top10us")
iv2 <- run_iv(panel, "ratio_5010")
run_heterogeneity(panel, "top10us")

p_plot <- plot_tradeoff(panel)
ggsave(paste0(DATA_DIR, "tradeoff_frontier.png"), p_plot,
       width=10, height=5.5, dpi=200)

cat("\nDone.\n")
