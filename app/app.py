import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
from pathlib import Path

PROC = Path("data/processed")
MIN_N = 30

st.set_page_config(page_title="AEI WardWatch", layout="wide")


@st.cache_data
def load():
    ab = pd.read_parquet(PROC / "antibiogram.parquet")
    sig = pd.read_parquet(PROC / "signals.parquet")
    wk = pd.read_parquet(PROC / "weekly.parquet")
    return ab, sig, wk


ab, sig, wk = load()

st.title("AEI WardWatch")
st.caption("Rolling, ward-level antibiotic efficacy indicator · last 4 weeks vs. annual baseline")

# ---------------------------------------------------------------- sidebar
with st.sidebar:
    st.header("Filters")
    wards = sorted(ab["ward"].astype(str).unique())
    default_ward_idx = wards.index("SURG") if "SURG" in wards else 0
    ward = st.selectbox("Ward", wards, index=default_ward_idx)

    specs = sorted(ab[ab["ward"].astype(str) == ward]["specimen"].astype(str).unique())
    target_specimen = "wound"
    default_spec_idx = specs.index(target_specimen) if target_specimen in specs else 0
    specimen = st.selectbox("Specimen", specs, index=default_spec_idx)

    orgs = sorted(ab[
        (ab["ward"].astype(str) == ward) & (ab["specimen"].astype(str) == specimen)
    ]["organism"].astype(str).unique())
    default_org_idx = orgs.index("Pseudomonas aeruginosa") if "Pseudomonas aeruginosa" in orgs else 0
    organism = st.selectbox("Organism", orgs, index=default_org_idx)
    period = st.radio("Period", ["4w", "annual"], horizontal=True,
                      format_func=lambda p: "Last 4 weeks" if p == "4w" else "Annual")

mask = (
    (ab["ward"].astype(str) == ward)
    & (ab["specimen"].astype(str) == specimen)
    & (ab["organism"].astype(str) == organism)
)

# ------------------------------------------------------- anomaly banner
ward_sig = sig[sig["ward"].astype(str) == ward]
if len(ward_sig):
    st.error(f"⚠ {len(ward_sig)} resistance signal(s) on {ward} in the last 4 weeks")
    with st.expander("Show signals"):
        st.dataframe(
            ward_sig[["organism", "antibiotic", "specimen", "n_4w", "pct_r_4w",
                      "n_base", "pct_r_base", "delta_pct", "z"]]
            .rename(columns={"pct_r_4w": "%R last 4w", "pct_r_base": "%R baseline",
                             "delta_pct": "Δ pp", "n_4w": "n 4w", "n_base": "n baseline"}),
            use_container_width=True, hide_index=True,
        )
else:
    st.success(f"No resistance signals on {ward} in the last 4 weeks")

# -------------------------------------------------------- antibiogram table
st.subheader(f"Antibiogram — {ward} · {specimen} · {organism}")

view = ab[mask & (ab["period"] == period)].copy()
if view.empty:
    st.warning(f"⚠ No {period} data available for this specific ward/organism combination. Try switching the Period to 'Annual'.")
else:
    # Vectorized string creation avoids pandas apply issues for empty selections.
    view["display"] = np.where(
        view["n"] >= MIN_N,
        view["pct_r"].astype(str) + "%  (n=" + view["n"].astype(int).astype(str) + ")",
        "Insufficient data (n=" + view["n"].astype(int).astype(str) + ")",
    )
    table = view[["antibiotic", "display"]].sort_values("antibiotic")
    st.dataframe(table.rename(columns={"display": "Resistance (%R)"}),
                 use_container_width=True, hide_index=True)
st.caption(f"Minimum {MIN_N} first isolates required for a reliable percentage.")

# --------------------------------------------------------------- trend
st.subheader("Weekly resistance trend")
ts = wk[(wk["ward"].astype(str) == ward) & (wk["specimen"].astype(str) == specimen)
        & (wk["organism"].astype(str) == organism)]
if len(ts):
    abx = st.selectbox("Antibiotic", sorted(ts["antibiotic"].astype(str).unique()))
    t = ts[ts["antibiotic"].astype(str) == abx].sort_values("week")
    fig = px.line(t, x="week", y="pct_r", markers=True,
                  labels={"week": "Week", "pct_r": "% resistant"})
    fig.update_layout(yaxis_range=[0, 100], height=320, margin=dict(t=10, b=10))
    st.plotly_chart(fig, use_container_width=True)
else:
    st.info("No weekly series for this selection.")