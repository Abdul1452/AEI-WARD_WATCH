"""
AEI WardWatch - data build.

1. Discovers raw CSVs and prints their schema (so we can verify assumptions).
2. Joins patients -> cultures -> isolates -> susceptibility.
3. Applies first-isolate-per-patient filtering.
4. Writes data/processed/{antibiogram,weekly,signals}.parquet
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd

RAW = os.path.join("data", "raw")
OUT = os.path.join("data", "processed")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- discovery
def find_csvs() -> list[str]:
    paths = sorted(glob.glob(os.path.join(RAW, "**", "*.csv"), recursive=True))
    if not paths:
        sys.exit(f"[FATAL] No CSVs under {RAW}/. Extract synthetic_microbiology_dataset.zip there.")
    return paths


def pick(cols, *patterns):
    """First column whose lowercased name contains any pattern."""
    low = {c: str(c).lower() for c in cols}
    for pat in patterns:
        for c, lc in low.items():
            if pat in lc:
                return c
    return None


def main() -> None:
    paths = find_csvs()
    print("\n=== RAW FILES ===")
    tables: dict[str, pd.DataFrame] = {}
    for p in paths:
        df = pd.read_csv(p, low_memory=False)
        name = os.path.splitext(os.path.basename(p))[0].lower()
        tables[name] = df
        print(f"  {name:24s} {len(df):>7,} rows x {df.shape[1]:>2} cols   {p}")

    print("\n=== COLUMNS ===")
    for name, df in tables.items():
        print(f"  {name}: {list(df.columns)}")

    # ------------------------------------------------------ classify tables
    # Find the susceptibility table: has an isolate link AND a result-like column,
    # or is in wide format (many columns whose values are S/I/R).
    def looks_like_sir(series: pd.Series) -> bool:
        vals = series.dropna().astype(str).str.upper().str[:1]
        return len(vals) > 0 and vals.isin(["S", "I", "R"]).mean() > 0.9

    sus_name = None
    for name, df in tables.items():
        if "isolate" not in " ".join(df.columns).lower():
            continue
        wide = [c for c in df.columns if looks_like_sir(df[c])]
        has_long = pick(df.columns, "antibiotic", "antimicrobial", "drug", "agent") is not None
        if has_long or len(wide) >= 3:
            sus_name = name
            break
    if sus_name is None:
        sys.exit("[FATAL] Could not identify the susceptibility table. Paste the COLUMNS output back.")

    iso_name = next((n for n, d in tables.items() if n != sus_name and "isolate_id" in d.columns), None)
    cul_name = next((n for n, d in tables.items()
                     if n not in (sus_name, iso_name) and "culture_id" in d.columns), None)
    if iso_name is None or cul_name is None:
        sys.exit("[FATAL] Could not identify isolates/cultures tables. Paste the COLUMNS output back.")

    print(f"\n  -> susceptibility: {sus_name}\n  -> isolates: {iso_name}\n  -> cultures: {cul_name}")

    sus = tables[sus_name].copy()
    iso = tables[iso_name].copy()
    cul = tables[cul_name].copy()

    # ------------------------------------------- normalise susceptibility
    ab_col = pick(sus.columns, "antibiotic", "antimicrobial", "drug", "agent")
    res_col = pick(sus.columns, "result", "susceptib", "interpret", "sir", "phenotype")

    if ab_col and res_col:
        sus = sus[[c for c in ["isolate_id", ab_col, res_col] if c in sus.columns]].rename(
            columns={ab_col: "antibiotic", res_col: "result"}
        )
    else:
        # wide format -> melt
        id_cols = [c for c in sus.columns if looks_like_sir(sus[c]) is False]
        value_cols = [c for c in sus.columns if looks_like_sir(sus[c])]
        print(f"  wide susceptibility detected, melting {len(value_cols)} antibiotic columns")
        sus = sus.melt(id_vars=id_cols, value_vars=value_cols,
                       var_name="antibiotic", value_name="result")

    # ------------------------------------------- normalise cultures/isolates
    date_col = pick(cul.columns, "date", "collect", "taken", "sample_d")
    ward_col = pick(cul.columns, "ward", "unit", "location", "clinic", "centre", "center")
    spec_col = pick(cul.columns, "specimen", "sample_type", "material", "sample")
    org_col = pick(iso.columns, "organism", "species", "bacteria", "pathogen", "microbe")

    missing = [n for n, v in [("date", date_col), ("ward", ward_col),
                              ("specimen", spec_col), ("organism", org_col)] if v is None]
    if missing:
        sys.exit(f"[FATAL] missing columns: {missing}. Paste the COLUMNS output back.")

    cul_small = cul[[c for c in ["culture_id", "patient_id", date_col, ward_col, spec_col] if c in cul.columns]]
    iso_small = iso[[c for c in ["isolate_id", "culture_id", org_col] if c in iso.columns]]

    df = (
        iso_small
        .merge(cul_small, on="culture_id", how="inner")
        .merge(sus, on="isolate_id", how="inner")
        .rename(columns={date_col: "date", ward_col: "ward",
                         spec_col: "specimen", org_col: "organism"})
    )
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df["result"] = df["result"].astype(str).str.upper().str.strip().str[:1]
    df = df[df["result"].isin(["S", "I", "R"])]
    df = df.dropna(subset=["date", "ward", "specimen", "organism", "antibiotic", "patient_id"])
    print(f"\n=== JOINED === {len(df):,} susceptibility records")

    # ------------------------------------- first isolate per patient+organism
    df = df.sort_values(["patient_id", "organism", "date"])
    keep_idx, last_seen = [], {}
    EPISODE_DAYS = 365
    for idx, pid, org, d in zip(df.index, df["patient_id"], df["organism"], df["date"]):
        key = (pid, org)
        prev = last_seen.get(key)
        if prev is None or (d - prev).days >= EPISODE_DAYS:
            keep_idx.append(idx)
            last_seen[key] = d
    df = df.loc[keep_idx]
    print(f"=== FIRST-ISOLATE FILTER === {len(df):,} records retained")

    # ------------------------------------------------------------ aggregate
    keys = ["ward", "specimen", "organism", "antibiotic"]
    max_date = df["date"].max()
    w4_start = max_date - pd.Timedelta(28, unit="D")
    y1_start = max_date - pd.Timedelta(365, unit="D")

    def agg(sub: pd.DataFrame, label: str) -> pd.DataFrame:
        g = sub.groupby(keys, dropna=False).agg(
            n=("result", "size"), n_r=("result", lambda s: int((s == "R").sum()))
        ).reset_index()
        g["pct_r"] = (g["n_r"] / g["n"] * 100).round(1)
        g["period"] = label
        return g

    last4w = agg(df[df["date"] > w4_start], "4w")
    annual = agg(df[df["date"] > y1_start], "annual")
    baseline = agg(df[(df["date"] > y1_start) & (df["date"] <= w4_start)], "baseline")

    antibiogram = pd.concat([last4w, annual], ignore_index=True)

    # weekly series (for trend chart)
    wk = df[df["date"] > y1_start].copy()
    wk["week"] = wk["date"].dt.to_period("W").dt.start_time
    weekly = (wk.groupby(keys + ["week"], dropna=False)
                .agg(n=("result", "size"), n_r=("result", lambda s: int((s == "R").sum())))
                .reset_index())
    weekly["pct_r"] = (weekly["n_r"] / weekly["n"] * 100).round(1)

    # ------------------------------------------------------------- signals
    sig = last4w.merge(baseline, on=keys, suffixes=("_4w", "_base"))
    
    # Relaxed thresholds to catch earlier signals
    sig = sig[(sig["n_4w"] >= 5) & (sig["n_base"] >= 20)].copy()
    
    p1 = sig["n_r_4w"] / sig["n_4w"]
    p2 = sig["n_r_base"] / sig["n_base"]
    p = (sig["n_r_4w"] + sig["n_r_base"]) / (sig["n_4w"] + sig["n_base"])
    se = np.sqrt(p * (1 - p) * (1 / sig["n_4w"] + 1 / sig["n_base"]))
    sig["z"] = ((p1 - p2) / se.replace(0, np.nan)).round(2)
    sig["delta_pct"] = (sig["pct_r_4w"] - sig["pct_r_base"]).round(1)
    
    # Show top 5 deviations in terminal so we can pick a good demo ward
    print("\n=== TOP 5 DEVIATIONS ===")
    print(sig[["ward", "organism", "antibiotic", "pct_r_4w", "pct_r_base", "delta_pct", "z"]]
          .sort_values("z", key=abs, ascending=False).head(5).to_string(index=False))
    
    # Filter for display: 90% confidence (|z| >= 1.64)
    sig = sig[sig["z"].abs() >= 1.64].sort_values("z", key=np.abs, ascending=False)
    print(f"=== SIGNALS === {len(sig)} strata flagged")

    # -------------------------------------------------------------- write
    antibiogram.to_parquet(os.path.join(OUT, "antibiogram.parquet"), index=False)
    weekly.to_parquet(os.path.join(OUT, "weekly.parquet"), index=False)
    sig.to_parquet(os.path.join(OUT, "signals.parquet"), index=False)
    print(f"\nWrote {len(antibiogram):,} antibiogram rows, {len(weekly):,} weekly rows, "
          f"{len(sig):,} signals to {OUT}/")
    print(f"Date range after filtering: {df['date'].min().date()} -> {max_date.date()}")
    print(f"Wards: {sorted(df['ward'].astype(str).unique())}")
    print(f"Organisms ({df['organism'].nunique()}): {sorted(df['organism'].astype(str).unique())[:10]} ...")


if __name__ == "__main__":
    main()