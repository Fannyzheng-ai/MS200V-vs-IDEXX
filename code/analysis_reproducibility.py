#!/usr/bin/env python3
"""
Reproduce the main-text tables of the MS200V vs IDEXX Catalyst One method
comparison study. Pearson r, Deming regression (variance ratio = 1, jackknife
95% CIs), and Bland-Altman percent-bias summary per analyte.

The CSV `data/paired_results.csv` is already filtered to the per-analyte
cohorts used in the published analysis; this script needs no further
exclusion logic.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import t as tdist

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "paired_results.csv")
TEA  = {"ALB": 15, "ALP": 25, "AMY": 25, "Ca": 10, "CHOL": 20, "CREA": 20,
        "GLU": 20, "PHOS": 15, "TBIL": 25, "TP": 10, "UREA": 15}
HYBRID_ABS = {"TBIL": 3.0}  # TBIL hybrid acceptance criterion


def deming(x, y, lam=1.0):
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x); mx, my = x.mean(), y.mean()
    sxx = np.sum((x-mx)**2)/(n-1); syy = np.sum((y-my)**2)/(n-1)
    sxy = np.sum((x-mx)*(y-my))/(n-1)
    A = syy - lam*sxx
    slope = (A + np.sqrt(A*A + 4*lam*sxy*sxy)) / (2*sxy)
    intercept = my - slope*mx
    r = sxy / np.sqrt(sxx*syy)
    return slope, intercept, r, n


def jackknife(x, y, alpha=0.05):
    n = len(x)
    slopes, ints = [], []
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        s, b, *_ = deming(x[m], y[m])
        slopes.append(s); ints.append(b)
    s_full, b_full, *_ = deming(x, y)
    slopes, ints = np.array(slopes), np.array(ints)
    se_s = np.sqrt((n-1)/n * np.sum((slopes - slopes.mean())**2))
    se_b = np.sqrt((n-1)/n * np.sum((ints - ints.mean())**2))
    tval = tdist.ppf(1 - alpha/2, n - 2)
    return s_full, (s_full - tval*se_s, s_full + tval*se_s), b_full, (b_full - tval*se_b, b_full + tval*se_b)


def main():
    df = pd.read_csv(DATA).dropna(subset=["IDEXX", "MS200V"])
    rows = []
    for a in sorted(df["Analyte"].unique()):
        d = df[df["Analyte"] == a]
        x = np.array(d["IDEXX"], float); y = np.array(d["MS200V"], float)
        slope, slope_ci, intercept, intercept_ci = jackknife(x, y)
        _, _, r, n = deming(x, y)
        pct = (y - x) / x * 100
        mean_bias, sd_bias = pct.mean(), pct.std(ddof=1)
        loa_lo, loa_hi = mean_bias - 1.96*sd_bias, mean_bias + 1.96*sd_bias
        tea = TEA[a]
        within = np.abs(pct) <= tea
        if a in HYBRID_ABS:
            within = within | (np.abs(y - x) <= HYBRID_ABS[a])
        pct_within = within.mean() * 100
        rows.append({
            "Analyte": a, "N": n,
            "Deming_slope": round(slope, 3),
            "Deming_slope_lo": round(slope_ci[0], 3),
            "Deming_slope_hi": round(slope_ci[1], 3),
            "Deming_intercept": round(intercept, 2),
            "Deming_intercept_lo": round(intercept_ci[0], 2),
            "Deming_intercept_hi": round(intercept_ci[1], 2),
            "Pearson_r": round(r, 3),
            "R2": round(r*r, 3),
            "Mean_bias_pct": round(mean_bias, 2),
            "SD_bias_pct": round(sd_bias, 2),
            "LoA_lo_pct": round(loa_lo, 1),
            "LoA_hi_pct": round(loa_hi, 1),
            "Pct_within_TEa": round(pct_within, 1),
            "ASVCP_TEa_pct": tea,
        })
    out = pd.DataFrame(rows)
    out_path = os.path.join(ROOT, "data", "expected_outputs.csv")
    out.to_csv(out_path, index=False)
    print(out.to_string(index=False))
    print(f"\nWritten: {out_path}")


if __name__ == "__main__":
    main()
