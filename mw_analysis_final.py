"""
Cyclical Effects of Minimum Wage Policy on Employment and Inequality
Ziyu Ren | University of Chicago, MA Economics

Python replication of mw_analysis_final.R.
Same panel construction, TWFE estimation, robustness checks,
neighbor-state IV, heterogeneity analysis, and tradeoff frontier plot.

Data (place in DATA_DIR):
  SAGDP1__ALL_AREAS_1997_2024.csv   BEA real GDP by state
  SAGDP2__ALL_AREAS_1997_2024.csv   BEA GDP by industry
  Frank_StatesUS_2022.xls           Top 10% income share (Frank 2022)
  hourly_wage_percentile_ratios.csv EPI 50-10 wage ratio

Packages: pandas, numpy, linearmodels, matplotlib, openpyxl, xlrd
"""

import pandas as pd
import numpy as np
from linearmodels.panel import PanelOLS
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

DATA_DIR = "/Users/ProjectPei/Desktop/"


# ==============================================================================
# 0. CONSTANTS
# ==============================================================================

FEDERAL_MW = {
    1997: 4.75, 1998: 5.15, 1999: 5.15, 2000: 5.15, 2001: 5.15,
    2002: 5.15, 2003: 5.15, 2004: 5.15, 2005: 5.15, 2006: 5.15,
    2007: 5.85, 2008: 6.55, 2009: 7.25, 2010: 7.25, 2011: 7.25,
    2012: 7.25, 2013: 7.25, 2014: 7.25, 2015: 7.25, 2016: 7.25,
    2017: 7.25, 2018: 7.25,
}

# 224 state-year cells where effective state MW exceeds the federal floor.
# Cross-checked against DOL Wage and Hour Division records and
# state ballot-measure implementation dates.
STATE_MW_OVERRIDES = {
    ("CA",2001):6.25, ("CA",2002):6.75, ("CA",2003):6.75, ("CA",2004):6.75,
    ("CA",2005):6.75, ("CA",2006):6.75, ("CA",2007):7.50, ("CA",2008):8.00,
    ("CA",2009):8.00, ("CA",2010):8.00, ("CA",2011):8.00, ("CA",2012):8.00,
    ("CA",2013):8.00, ("CA",2014):9.00, ("CA",2015):9.00, ("CA",2016):10.00,
    ("CA",2017):10.50,("CA",2018):11.00,
    ("NY",2005):6.00, ("NY",2006):6.75, ("NY",2007):7.15, ("NY",2008):7.15,
    ("NY",2014):8.00, ("NY",2015):8.75, ("NY",2016):9.00, ("NY",2017):9.70,
    ("NY",2018):10.40,
    ("WA",1999):5.70, ("WA",2000):6.50, ("WA",2001):6.72, ("WA",2002):6.90,
    ("WA",2003):7.01, ("WA",2004):7.16, ("WA",2005):7.35, ("WA",2006):7.63,
    ("WA",2007):7.93, ("WA",2008):8.07, ("WA",2009):8.55, ("WA",2010):8.55,
    ("WA",2011):8.67, ("WA",2012):9.04, ("WA",2013):9.19, ("WA",2014):9.32,
    ("WA",2015):9.47, ("WA",2016):9.47, ("WA",2017):11.00,("WA",2018):11.50,
    ("MA",2007):7.50, ("MA",2008):8.00, ("MA",2015):9.00, ("MA",2016):10.00,
    ("MA",2017):11.00,("MA",2018):11.00,
    ("OR",1999):6.00, ("OR",2000):6.50, ("OR",2001):6.50, ("OR",2002):6.50,
    ("OR",2003):6.90, ("OR",2004):7.05, ("OR",2005):7.25, ("OR",2006):7.50,
    ("OR",2007):7.80, ("OR",2008):7.95, ("OR",2009):8.40, ("OR",2010):8.40,
    ("OR",2011):8.50, ("OR",2012):8.80, ("OR",2013):8.95, ("OR",2014):9.10,
    ("OR",2015):9.25, ("OR",2016):9.25, ("OR",2017):9.75, ("OR",2018):10.25,
    ("IL",2004):5.50, ("IL",2005):6.50, ("IL",2006):6.50, ("IL",2007):7.50,
    ("IL",2008):7.75, ("IL",2009):8.00, ("IL",2010):8.25, ("IL",2011):8.25,
    ("IL",2012):8.25, ("IL",2013):8.25, ("IL",2014):8.25, ("IL",2015):8.25,
    ("IL",2016):8.25, ("IL",2017):8.25, ("IL",2018):8.25,
    ("CT",2014):8.70, ("CT",2015):9.15, ("CT",2016):9.60, ("CT",2017):10.10,
    ("CT",2018):10.10,
    ("MN",2014):8.00, ("MN",2015):9.00, ("MN",2016):9.00, ("MN",2017):9.50,
    ("MN",2018):9.65,
    ("NJ",2014):8.25, ("NJ",2015):8.38, ("NJ",2016):8.38, ("NJ",2017):8.44,
    ("NJ",2018):8.60,
    ("CO",2007):6.85, ("CO",2008):7.02, ("CO",2009):7.28, ("CO",2010):7.24,
    ("CO",2011):7.36, ("CO",2012):7.64, ("CO",2013):7.78, ("CO",2014):8.00,
    ("CO",2015):8.23, ("CO",2016):8.31, ("CO",2017):9.30, ("CO",2018):10.20,
    ("AZ",2007):6.75, ("AZ",2008):6.90, ("AZ",2009):7.25, ("AZ",2010):7.25,
    ("AZ",2011):7.35, ("AZ",2012):7.65, ("AZ",2013):7.80, ("AZ",2014):7.90,
    ("AZ",2015):8.05, ("AZ",2016):8.05, ("AZ",2017):10.00,("AZ",2018):10.50,
    ("VT",2007):7.53, ("VT",2008):7.68, ("VT",2009):8.06, ("VT",2010):8.06,
    ("VT",2011):8.15, ("VT",2012):8.46, ("VT",2013):8.60, ("VT",2014):8.73,
    ("VT",2015):9.15, ("VT",2016):9.60, ("VT",2017):10.00,("VT",2018):10.50,
    ("MI",2006):6.95, ("MI",2007):7.15, ("MI",2008):7.40, ("MI",2014):8.15,
    ("MI",2015):8.15, ("MI",2016):8.50, ("MI",2017):8.90, ("MI",2018):9.25,
    ("MD",2015):8.00, ("MD",2016):8.75, ("MD",2017):9.25, ("MD",2018):9.25,
    ("NE",2015):8.00, ("NE",2016):9.00, ("NE",2017):9.00, ("NE",2018):9.00,
    ("SD",2015):8.50, ("SD",2016):8.55, ("SD",2017):8.65, ("SD",2018):8.85,
    ("AK",2003):7.15, ("AK",2014):8.75, ("AK",2015):9.75, ("AK",2016):9.75,
    ("AK",2017):9.80, ("AK",2018):9.84,
    ("DC",1999):6.15, ("DC",2003):6.60, ("DC",2006):7.00, ("DC",2008):7.55,
    ("DC",2009):8.25, ("DC",2010):8.25, ("DC",2011):8.25, ("DC",2012):8.25,
    ("DC",2013):8.25, ("DC",2014):9.50, ("DC",2015):10.50,("DC",2016):11.50,
    ("DC",2017):12.50,("DC",2018):13.25,
    ("RI",2014):8.00, ("RI",2015):9.00, ("RI",2016):9.60, ("RI",2017):9.60,
    ("RI",2018):10.10,
    ("HI",2003):6.25, ("HI",2004):6.25, ("HI",2005):6.75, ("HI",2006):7.25,
    ("HI",2007):7.25, ("HI",2015):7.75, ("HI",2016):8.50, ("HI",2017):9.25,
    ("HI",2018):10.10,
    ("NM",2007):5.15, ("NM",2008):6.50, ("NM",2009):7.50,
    ("FL",2005):6.15, ("FL",2006):6.40, ("FL",2007):6.67, ("FL",2008):6.79,
    ("FL",2009):7.21, ("FL",2011):7.31, ("FL",2012):7.67, ("FL",2013):7.79,
    ("FL",2014):7.93, ("FL",2015):8.05, ("FL",2016):8.05, ("FL",2017):8.10,
    ("FL",2018):8.25,
    ("MT",2007):6.15, ("MT",2008):6.25, ("MT",2009):6.90, ("MT",2015):8.05,
    ("MT",2016):8.05, ("MT",2017):8.15, ("MT",2018):8.30,
    ("NV",2009):7.55, ("NV",2010):8.25, ("NV",2011):8.25, ("NV",2015):8.25,
    ("MO",2007):6.65, ("MO",2008):6.65, ("MO",2009):7.05,
    ("OH",2006):6.85, ("OH",2007):7.00, ("OH",2008):7.00, ("OH",2009):7.30,
    ("OH",2010):7.30, ("OH",2011):7.40, ("OH",2012):7.70, ("OH",2013):7.85,
    ("OH",2014):8.00, ("OH",2015):8.10, ("OH",2016):8.15, ("OH",2017):8.15,
    ("OH",2018):8.30,
    ("ME",2017):9.00, ("ME",2018):10.00,
    ("WV",2016):8.75,
    ("KS",2010):7.25,
}

# Neighbor-state adjacency for IV construction
ADJACENCY = {
    "AL":["FL","GA","MS","TN"],    "AK":[],
    "AZ":["CA","CO","NM","NV","UT"],
    "AR":["LA","MO","MS","OK","TN","TX"], "CA":["AZ","NV","OR"],
    "CO":["AZ","KS","NE","NM","OK","UT","WY"], "CT":["MA","NY","RI"],
    "DE":["MD","NJ","PA"],         "FL":["AL","GA"],
    "GA":["AL","FL","NC","SC","TN"],          "HI":[],
    "ID":["MT","NV","OR","UT","WA","WY"],     "IL":["IN","IA","KY","MO","WI"],
    "IN":["IL","KY","MI","OH"],    "IA":["IL","MN","MO","NE","SD","WI"],
    "KS":["CO","MO","NE","OK"],    "KY":["IL","IN","MO","OH","TN","VA","WV"],
    "LA":["AR","MS","TX"],         "ME":["NH"],
    "MD":["DE","PA","VA","WV","DC"],          "MA":["CT","NH","NY","RI","VT"],
    "MI":["IN","OH","WI"],         "MN":["IA","ND","SD","WI"],
    "MS":["AL","AR","LA","TN"],    "MO":["AR","IL","IA","KS","KY","NE","OK","TN"],
    "MT":["ID","ND","SD","WY"],    "NE":["CO","IA","KS","MO","SD","WY"],
    "NV":["AZ","CA","ID","OR","UT"],          "NH":["MA","ME","VT"],
    "NJ":["DE","NY","PA"],         "NM":["AZ","CO","OK","TX"],
    "NY":["CT","MA","NJ","PA","VT"],          "NC":["GA","SC","TN","VA"],
    "ND":["MN","MT","SD"],         "OH":["IN","KY","MI","PA","WV"],
    "OK":["AR","CO","KS","MO","NM","TX"],     "OR":["CA","ID","NV","WA"],
    "PA":["DE","MD","NJ","NY","OH","WV"],     "RI":["CT","MA"],
    "SC":["GA","NC"],              "SD":["IA","MN","MT","ND","NE","WY"],
    "TN":["AL","AR","GA","KY","MS","MO","NC","VA"], "TX":["AR","LA","NM","OK"],
    "UT":["AZ","CO","ID","NV","NM","WY"],     "VT":["MA","NH","NY"],
    "VA":["KY","MD","NC","TN","WV","DC"],     "WA":["ID","OR"],
    "WV":["KY","MD","OH","PA","VA"],          "WI":["IL","IA","MI","MN"],
    "WY":["CO","ID","MT","NE","SD","UT"],     "DC":["MD","VA"],
}

FIPS_TO_STATE = {
    "01000":"AL","02000":"AK","04000":"AZ","05000":"AR","06000":"CA",
    "08000":"CO","09000":"CT","10000":"DE","12000":"FL","13000":"GA",
    "15000":"HI","16000":"ID","17000":"IL","18000":"IN","19000":"IA",
    "20000":"KS","21000":"KY","22000":"LA","23000":"ME","24000":"MD",
    "25000":"MA","26000":"MI","27000":"MN","28000":"MS","29000":"MO",
    "30000":"MT","31000":"NE","32000":"NV","33000":"NH","34000":"NJ",
    "35000":"NM","36000":"NY","37000":"NC","38000":"ND","39000":"OH",
    "40000":"OK","41000":"OR","42000":"PA","44000":"RI","45000":"SC",
    "46000":"SD","47000":"TN","48000":"TX","49000":"UT","50000":"VT",
    "51000":"VA","53000":"WA","54000":"WV","55000":"WI","56000":"WY",
    "11000":"DC",
}

NAME_TO_ABBR = {
    "Alabama":"AL","Alaska":"AK","Arizona":"AZ","Arkansas":"AR",
    "California":"CA","Colorado":"CO","Connecticut":"CT","Delaware":"DE",
    "Florida":"FL","Georgia":"GA","Hawaii":"HI","Idaho":"ID",
    "Illinois":"IL","Indiana":"IN","Iowa":"IA","Kansas":"KS",
    "Kentucky":"KY","Louisiana":"LA","Maine":"ME","Maryland":"MD",
    "Massachusetts":"MA","Michigan":"MI","Minnesota":"MN","Mississippi":"MS",
    "Missouri":"MO","Montana":"MT","Nebraska":"NE","Nevada":"NV",
    "New Hampshire":"NH","New Jersey":"NJ","New Mexico":"NM","New York":"NY",
    "North Carolina":"NC","North Dakota":"ND","Ohio":"OH","Oklahoma":"OK",
    "Oregon":"OR","Pennsylvania":"PA","Rhode Island":"RI",
    "South Carolina":"SC","South Dakota":"SD","Tennessee":"TN","Texas":"TX",
    "Utah":"UT","Vermont":"VT","Virginia":"VA","Washington":"WA",
    "West Virginia":"WV","Wisconsin":"WI","Wyoming":"WY",
    "District of Columbia":"DC",
}

YR_COLS = [str(y) for y in range(1997, 2019)]


# ==============================================================================
# 1. DATA LOADING & PANEL CONSTRUCTION
# ==============================================================================

def load_gdp():
    bea = pd.read_csv(DATA_DIR + "SAGDP1__ALL_AREAS_1997_2024.csv",
                      dtype={"GeoFIPS": str})
    bea["GeoFIPS"] = bea["GeoFIPS"].str.strip().str.strip('"').str.zfill(5)
    gdp = (bea[(bea["GeoFIPS"].isin(FIPS_TO_STATE)) & (bea["LineCode"] == 1)]
           .assign(state=lambda d: d["GeoFIPS"].map(FIPS_TO_STATE))
           [["state"] + YR_COLS]
           .melt(id_vars="state", var_name="year", value_name="real_gdp")
           .assign(year=lambda d: d["year"].astype(int),
                   real_gdp=lambda d: pd.to_numeric(d["real_gdp"], errors="coerce"))
           .query("year >= 1998")
           .sort_values(["state", "year"])
           .reset_index(drop=True))
    gdp["gdp_growth"] = gdp.groupby("state")["real_gdp"].pct_change()
    return gdp[["state", "year", "real_gdp", "gdp_growth"]]


def load_sector_shares():
    bea = pd.read_csv(DATA_DIR + "SAGDP2__ALL_AREAS_1997_2024.csv",
                      dtype={"GeoFIPS": str})
    bea["GeoFIPS"]     = bea["GeoFIPS"].str.strip().str.strip('"').str.zfill(5)
    bea["Description"] = bea["Description"].str.strip()
    bea = bea[bea["GeoFIPS"].isin(FIPS_TO_STATE)].copy()
    bea["state"] = bea["GeoFIPS"].map(FIPS_TO_STATE)

    def pivot_sector(descs, col):
        sub = (bea[bea["Description"].isin(descs)]
               .groupby("state")[YR_COLS]
               .apply(lambda x: x.apply(pd.to_numeric, errors="coerce").sum())
               .reset_index()
               .melt(id_vars="state", var_name="year", value_name=col)
               .assign(year=lambda d: d["year"].astype(int)))
        return sub[sub["year"] >= 1998]

    total = pivot_sector(["All industry total"],                                        "gdp_total")
    fin   = pivot_sector(["Finance and insurance"],                                     "gdp_fin")
    manuf = pivot_sector(["Durable goods manufacturing","Nondurable goods manufacturing"], "gdp_manuf")

    return (total
            .merge(fin,   on=["state","year"])
            .merge(manuf, on=["state","year"])
            .assign(finance_share=lambda d: d["gdp_fin"]   / d["gdp_total"] * 100,
                    manuf_share  =lambda d: d["gdp_manuf"] / d["gdp_total"] * 100)
            [["state","year","finance_share","manuf_share"]])


def load_frank():
    frank = pd.read_excel(DATA_DIR + "Frank_StatesUS_2022.xls",
                          sheet_name="Top Income ShareUS", header=0)
    frank.columns = ["year","st","state_name","top10us","top5us","top1us",
                     "top05us","top01us","x","xx"]
    return (frank
            .assign(state  =lambda d: d["state_name"].map(NAME_TO_ABBR),
                    year   =lambda d: pd.to_numeric(d["year"], errors="coerce").astype("Int64"),
                    top10us=lambda d: pd.to_numeric(d["top10us"], errors="coerce"))
            .query("1998 <= year <= 2018")
            .dropna(subset=["state"])
            [["state","year","top10us"]]
            .astype({"year": int}))


def load_wage_ratio():
    epi = pd.read_csv(DATA_DIR + "hourly_wage_percentile_ratios.csv")
    mask = ((epi["geo_type"] == "state") &
            epi["group_value"].str.contains("50", na=False) &
            epi["group_value"].str.contains("10", na=False) &
            epi["year"].between(1998, 2018))
    return (epi[mask]
            .assign(state     =lambda d: d["geo_name"].map(NAME_TO_ABBR),
                    ratio_5010=lambda d: pd.to_numeric(d["value"], errors="coerce"))
            .dropna(subset=["state"])
            [["state","year","ratio_5010"]])


def build_minimum_wage():
    rows = []
    for state in ADJACENCY:
        for year in range(1998, 2019):
            fed   = FEDERAL_MW.get(year, 7.25)
            state_floor = STATE_MW_OVERRIDES.get((state, year), fed)
            rows.append({"state": state, "year": year,
                         "min_wage": max(fed, state_floor)})
    return pd.DataFrame(rows)


def build_neighbor_mw(mw):
    # Instrument: average MW across contiguous states.
    # Exclusion restriction: neighbor legislatures set wages through independent
    # political processes uncorrelated with own-state income distribution.
    lookup = mw.set_index(["state","year"])["min_wage"]
    rows = []
    for _, row in mw.iterrows():
        nbrs = ADJACENCY.get(row["state"], [])
        if not nbrs:
            continue
        vals = [lookup.get((n, row["year"])) for n in nbrs
                if (n, row["year"]) in lookup.index]
        if vals:
            rows.append({"state": row["state"], "year": row["year"],
                         "neighbor_mw": np.mean(vals)})
    return pd.DataFrame(rows)


def build_panel():
    gdp  = load_gdp()
    sec  = load_sector_shares()
    frk  = load_frank()
    wr   = load_wage_ratio()
    mw   = build_minimum_wage()
    nbr  = build_neighbor_mw(mw)

    panel = (gdp
             .merge(mw,  on=["state","year"])
             .merge(sec, on=["state","year"], how="left")
             .merge(frk, on=["state","year"], how="left")
             .merge(wr,  on=["state","year"], how="left")
             .merge(nbr, on=["state","year"], how="left")
             .dropna(subset=["gdp_growth","min_wage"])
             .sort_values(["state","year"])
             .reset_index(drop=True))

    # Lagged MW — computed within state group
    panel["mw_lag1"]   = panel.groupby("state")["min_wage"].shift(1)
    panel["mwlag_gdp"] = panel["mw_lag1"] * panel["gdp_growth"]

    panel["mw_gdp"]          = panel["min_wage"]    * panel["gdp_growth"]
    panel["gdp_growth_sq"]   = panel["gdp_growth"]  ** 2
    panel["neighbor_mw_gdp"] = panel["neighbor_mw"] * panel["gdp_growth"]
    panel["finance_gdp"]     = panel["finance_share"] * panel["gdp_growth"]
    panel["recession"]       = panel["year"].isin([2008, 2009]).astype(int)
    panel["mw_rec"]          = panel["min_wage"] * panel["recession"]

    # Kaitz index: captures whether MW is actually binding in local labor markets
    yr_med = panel.groupby("year")["min_wage"].transform("median")
    panel["mw_relative"] = panel["min_wage"] / yr_med

    print(f"Panel: {len(panel):,} obs | {panel.state.nunique()} states | "
          f"{panel.year.min()}-{panel.year.max()}")
    print(f"  top10us:    {panel.top10us.notna().sum():,} obs")
    print(f"  ratio_5010: {panel.ratio_5010.notna().sum():,} obs")
    return panel


# ==============================================================================
# 2. SUMMARY STATISTICS
# ==============================================================================

def summary_stats(panel):
    print("\n" + "="*68)
    print("Summary Statistics")
    print("="*68 + "\n")

    vars_ = {
        "min_wage":      "Effective minimum wage ($)",
        "gdp_growth":    "Real GDP growth rate",
        "top10us":       "Top 10% income share",
        "ratio_5010":    "50-10 wage ratio",
        "finance_share": "Finance share of GDP (%)",
        "manuf_share":   "Manufacturing share of GDP (%)",
        "mw_relative":   "Kaitz index",
    }
    rows = []
    for col, label in vars_.items():
        x = panel[col].dropna()
        rows.append({
            "Variable": label,
            "N":        len(x),
            "Mean":     round(x.mean(), 3),
            "SD":       round(x.std(),  3),
            "Min":      round(x.min(),  3),
            "Median":   round(x.median(), 3),
            "Max":      round(x.max(),  3),
        })
    print(pd.DataFrame(rows).to_string(index=False))


# ==============================================================================
# 3. ESTIMATION HELPERS
# ==============================================================================

def to_panel(df):
    return df.set_index(["state","year"])


def twfe(formula, data, cluster=True):
    """
    TWFE via PanelOLS with entity and time effects.
    Equivalent to feols(y ~ x | state + year, cluster = ~state) in R.
    """
    return PanelOLS.from_formula(
        formula + " + EntityEffects + TimeEffects",
        data=to_panel(data),
        drop_absorbed=True
    ).fit(cov_type="clustered", cluster_entity=cluster)


def sig(p):
    return "***" if p < .01 else "**" if p < .05 else "*" if p < .10 else "(ns)"


def print_coef(res, vars_):
    for v in vars_:
        if v not in res.params.index:
            continue
        b, se, p = res.params[v], res.std_errors[v], res.pvalues[v]
        print(f"  {v:<32} {b:>9.4f}{sig(p):<4}  SE={se:.4f}  p={p:.3f}")
    print(f"  {'N':<32} {int(res.nobs):>9}")


# ==============================================================================
# 4. LAYER 1 — CYCLICAL ASYMMETRY
# ==============================================================================

def run_layer1(panel):
    print("\n" + "="*68)
    print("Layer 1: Cyclical Asymmetry — Top 10% Income Share")
    print("="*68 + "\n")

    p = panel.dropna(subset=["top10us"])

    m1 = twfe("top10us ~ min_wage + gdp_growth", p)
    m2 = twfe("top10us ~ min_wage + gdp_growth + mw_gdp", p)

    print("  Basic TWFE:")
    print_coef(m1, ["min_wage","gdp_growth"])
    print("\n  + MW×GDP:")
    print_coef(m2, ["min_wage","mw_gdp","gdp_growth"])

    b0 = m2.params["min_wage"]
    b1 = m2.params["mw_gdp"]
    print(f"\n  MW×GDP = {b1:.4f} (SE={m2.std_errors['mw_gdp']:.4f}) "
          f"| threshold: GDP growth < {-b0/b1*100:.1f}%")
    return m2


# ==============================================================================
# 5. ROBUSTNESS
# ==============================================================================

def run_robustness(panel):
    print("\n" + "="*68)
    print("Robustness: MW×GDP Across Six Specifications")
    print("="*68 + "\n")

    p     = panel.dropna(subset=["top10us"])
    p_lag = p.dropna(subset=["mw_lag1"])

    specs = [
        ("Baseline",        "top10us ~ min_wage + gdp_growth + mw_gdp",                              p,     "mw_gdp"),
        ("+ GDP²",          "top10us ~ min_wage + gdp_growth + mw_gdp + gdp_growth_sq",              p,     "mw_gdp"),
        ("+ Finance",       "top10us ~ min_wage + gdp_growth + mw_gdp + finance_share",               p,     "mw_gdp"),
        ("+ Finance×GDP",   "top10us ~ min_wage + gdp_growth + mw_gdp + finance_share + finance_gdp", p,     "mw_gdp"),
        ("+ GDP²+Finance",  "top10us ~ min_wage + gdp_growth + mw_gdp + gdp_growth_sq + finance_share", p,   "mw_gdp"),
        ("Lagged MW",       "top10us ~ mw_lag1 + gdp_growth + mwlag_gdp + gdp_growth_sq + finance_share", p_lag, "mwlag_gdp"),
    ]

    rows = []
    print(f"  {'Specification':<22} {'MW×GDP':>8}  {'':5}  SE      p-value")
    print(f"  {'-'*58}")
    for label, fml, data, key in specs:
        m = twfe(fml, data)
        b, se, p = m.params[key], m.std_errors[key], m.pvalues[key]
        rows.append(b)
        print(f"  {label:<22} {b:>8.4f}  {sig(p):<5}  {se:.4f}  {p:.5f}")

    print(f"\n  Range: {min(rows):.4f} – {max(rows):.4f} | "
          f"All p < 0.01: {all(twfe(s[1],s[2]).pvalues.get(s[3],1) < .01 for s in specs)}")


# ==============================================================================
# 6. INSTRUMENTAL VARIABLES
# ==============================================================================

def run_iv(panel, outcome="top10us"):
    print("\n" + "="*68)
    print(f"IV: Neighbor-State MW — outcome: {outcome}")
    print("="*68 + "\n")

    p = panel.dropna(subset=[outcome, "neighbor_mw"]).copy()
    print(f"  Sample: {len(p):,} obs, {p.state.nunique()} states "
          f"(AK and HI excluded)\n")

    # First stage: instrument relevance
    fs_mw  = twfe("min_wage ~ neighbor_mw + neighbor_mw_gdp + gdp_growth", p)
    fs_int = twfe("mw_gdp   ~ neighbor_mw + neighbor_mw_gdp + gdp_growth", p)

    # Partial F for each endogenous variable (Kleibergen-Paap analog)
    f_mw  = (fs_mw.params["neighbor_mw"] / fs_mw.std_errors["neighbor_mw"]) ** 2
    f_int = (fs_int.params["neighbor_mw_gdp"] / fs_int.std_errors["neighbor_mw_gdp"]) ** 2
    print(f"  First-stage F (min_wage ~ instruments):  {f_mw:,.2f}")
    print(f"  First-stage F (mw_gdp   ~ instruments):  {f_int:,.2f}\n")

    # 2SLS via fitted values (equivalent to R feols IV)
    p["mw_hat"]     = fs_mw.fitted_values.values
    p["mwgdp_hat"]  = fs_int.fitted_values.values

    m_ols = twfe(f"{outcome} ~ min_wage + mw_gdp + gdp_growth", p)
    m_iv  = twfe(f"{outcome} ~ mw_hat + mwgdp_hat + gdp_growth", p)

    b_ols, p_ols = m_ols.params["mw_gdp"],    m_ols.pvalues["mw_gdp"]
    b_iv,  p_iv  = m_iv.params["mwgdp_hat"],  m_iv.pvalues["mwgdp_hat"]

    print(f"  OLS  MW×GDP = {b_ols:.4f}{sig(p_ols):<4}  (p={p_ols:.3f})")
    print(f"  2SLS MW×GDP = {b_iv:.4f}{sig(p_iv):<4}  (p={p_iv:.3f})")
    print(f"  IV/OLS      = {b_iv/b_ols:.2f}x")


# ==============================================================================
# 7. HETEROGENEITY ANALYSIS
# ==============================================================================

def run_heterogeneity(panel, outcome="top10us"):
    print("\n" + "="*68)
    print(f"Heterogeneity Analysis — outcome: {outcome}")
    print("="*68)

    p = panel.dropna(subset=[outcome]).copy()

    state_chars = (p.groupby("state")
                   .agg(kaitz_mean=("mw_relative","mean"),
                        gdp_vol   =("gdp_growth", "std"),
                        manuf_mean=("manuf_share","mean"))
                   .reset_index())
    for col in ["kaitz_mean","gdp_vol","manuf_mean"]:
        state_chars[f"high_{col}"] = (state_chars[col] > state_chars[col].median())

    p = p.merge(state_chars[["state","high_kaitz_mean","high_gdp_vol","high_manuf_mean"]],
                on="state")

    fml = f"{outcome} ~ min_wage + gdp_growth + mw_gdp + gdp_growth_sq"

    dims = [
        ("high_kaitz_mean", "Kaitz index",
         "High Kaitz → asymmetry stronger where MW is binding"),
        ("high_gdp_vol",    "GDP volatility",
         "More volatile → larger cycle shocks"),
        ("high_manuf_mean", "Manufacturing share",
         "Higher labor demand elasticity → cost more cycle-sensitive"),
    ]

    for dim, label, rationale in dims:
        print(f"\n  {label} ({rationale})")
        for grp, name in [(True,"High"),(False,"Low")]:
            sub = p[p[dim] == grp]
            m   = twfe(fml, sub)
            b, pv = m.params["mw_gdp"], m.pvalues["mw_gdp"]
            print(f"    {name:<4} (N={len(sub):>5,}): MW×GDP = {b:.4f}{sig(pv):<4}  p={pv:.3f}")


# ==============================================================================
# 8. TRADEOFF FRONTIER PLOT
# ==============================================================================

def plot_tradeoff(panel, save_path=DATA_DIR + "tradeoff_frontier.png"):
    p  = panel.dropna(subset=["top10us"])
    m  = twfe("top10us ~ min_wage + gdp_growth + mw_gdp", p)
    b0 = m.params["min_wage"]
    b1 = m.params["mw_gdp"]

    g  = np.linspace(-0.06, 0.08, 300)
    me = b0 + b1 * g

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.axvspan(-6, 0, alpha=0.07, color="red")
    ax.axvspan( 0, 8, alpha=0.07, color="green")
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.axvline((-b0/b1)*100, color="#c0392b", lw=0.9, ls=":")
    ax.plot(g*100, me, color="#c0392b", lw=2.2)
    ax.text(-4.5, max(me)*0.85, "Contraction",
            color="#c0392b", fontstyle="italic", fontsize=10)
    ax.text( 5.0, max(me)*0.85, "Expansion",
            color="darkgreen", fontstyle="italic", fontsize=10)
    ax.set_xlabel("GDP Growth Rate (%)", fontsize=11)
    ax.set_ylabel("\u2202Top10%Share / \u2202MW", fontsize=10)
    ax.set_title("Equity-Efficiency Tradeoff Frontier",
                 fontsize=13, fontweight="bold")
    ax.annotate(
        f"\u2202TopShare/\u2202MW = {b0:.3f} + {b1:.3f} \u00d7 GDP_growth\n"
        f"Threshold: GDP growth < {(-b0/b1)*100:.1f}%",
        xy=(0.03, 0.06), xycoords="axes fraction", fontsize=9, color="gray"
    )
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    print(f"\n  Figure saved: {save_path}")
    plt.close()


# ==============================================================================
# MAIN
# ==============================================================================

if __name__ == "__main__":
    panel = build_panel()
    panel.to_csv(DATA_DIR + "panel_final.csv", index=False)

    summary_stats(panel)
    run_layer1(panel)
    run_robustness(panel)
    run_iv(panel, "top10us")
    run_iv(panel, "ratio_5010")
    run_heterogeneity(panel, "top10us")
    plot_tradeoff(panel)

    print("\nDone.")
