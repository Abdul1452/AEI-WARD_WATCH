# AEI WardWatch

**A rolling, real-time, ward-level antibiotic efficacy indicator (AEI) for antimicrobial resistance surveillance.**

Built for the **Tyks Laboratories** challenge at **JunctionX Turku 2026**.

> *A doctor choosing an antibiotic for a patient today cannot see the last four weeks of resistance data for their own ward. Annual, hospital-level reports are too slow and too coarse. WardWatch closes that gap.*

---

## The Problem

The WHO, UN, and World Economic Forum identify antimicrobial resistance (AMR) as one of the most serious global public health threats. Tyks Laboratories performs antimicrobial susceptibility testing and monitors bacterial resistance — but reports are currently compiled **annually at hospital level**.

That cadence is too infrequent and too coarse for the moment that matters most: **a clinician selecting an antibiotic for a patient on a specific ward, today**.

## What WardWatch Does

A single-screen tool for clinicians that answers:

> *"For my ward, over the last four weeks, what is the resistance profile of the organism I am treating — and has anything shifted?"*

**Core views:**
- **Ward-level antibiogram** with `Last 4 weeks` / `Annual` toggle
- **Cascading filters:** Ward → Specimen type → Organism
- **Anomaly signal banner** that surfaces statistically defensible deviations from the ward's annual baseline
- **Suppression logic** so no percentage is shown without a sufficient isolate count

## Methodology — Why the Numbers Are Trustworthy

Clinical surveillance tools fail when the underlying statistics are naive. WardWatch applies two established AMR-surveillance rules:

1. **First-isolate-per-patient filtering.** Repeat isolates from the same patient inflate resistance estimates. We include only the first isolate per patient per species within an episode window, aligning with EARS-Net methodology.
2. **Minimum 30-isolate threshold.** Strata below 30 isolates are shown as *"Insufficient data (n=X)"* rather than as a misleading percentage.

**Anomaly detection:** For each ward–organism–antibiotic stratum, the rolling 4-week resistance proportion is compared against that stratum's annual baseline. A signal is raised when the deviation exceeds a calibrated threshold. The signal is always displayed alongside its `n` and baseline so the clinician can judge credibility.

## Architecture

```
CSV dataset (synthetic microbiology dataset)
    ↓
Python / Pandas preprocessing (offline, reproducible)
    ├── Join patients → cultures → isolates → susceptibility
    ├── Apply first-isolate-per-patient filter
    └── Pre-aggregate into parquet (4-week rolling + annual strata)
    ↓
Streamlit app (loads parquet in milliseconds)
    ├── Pillar 1: ward antibiogram with filters
    └── Pillar 2: anomaly signal banner + detail view
```

**Stack:** Python, Pandas, PyArrow (parquet), Streamlit, Plotly.

## Getting Started

```bash
# 1. Clone
git clone https://github.com/<your-username>/aei-wardwatch.git
cd aei-wardwatch

# 2. Environment
python -m venv .venv
source .venv/bin/activate        # Windows: .\.venv\Scripts\Activate.ps1
pip install -r requirements.txt

# 3. Extract dataset
# Place synthetic_microbiology_dataset.zip in data/raw/ and extract it there.

# 4. Preprocess (produces data/processed/antibiogram.parquet)
python scripts/preprocess.py

# 5. Run the app
streamlit run app/app.py
```

App opens at `http://localhost:8501`.

## Project Structure

```
aei-wardwatch/
├── data/
│   ├── raw/              # extracted dataset (not committed)
│   └── processed/        # generated parquet (not committed)
├── scripts/
│   ├── inspect_data.py   # dataset structure explorer
│   └── preprocess.py     # join + first-isolate filter + aggregation
├── app/
│   └── app.py            # Streamlit UI
├── requirements.txt
└── README.md
```

## Dataset

Synthetic microbiology dataset provided by the challenge: ~10,000 patients, ~47,000 cultures (urine, blood, wound), three years, nine wards/health centres. Structure: `patients → cultures → isolates → susceptibility`, joined on `patient_id`, `culture_id`, `isolate_id`.

## Team

- Abdulwasiu Abdullahi Olamilekan — Software Engineer

## Hackathon Context

- **Event:** JunctionX Turku 2026 (16–18 October 2026)
- **Challenge:** Antibiotic efficacy indicator (AEI) — Tyks Laboratories
- **Pitch video:** `<link>` *(added at submission)*
- **Demo:** `<link or "run locally via instructions above">`

## References

- EARS-Net (ECDC) — European antimicrobial resistance surveillance: https://www.ecdc.europa.eu/en/about-us/networks/disease-networks-and-laboratory-networks/ears-net-data
- THL — Finnish antimicrobial resistance data: https://thl.fi/en/data-and-statistics/infectious-diseases/antimicrobial-resistance

## License

MIT — see `LICENSE`. (Junction rules confirm participants retain ownership of work created at the event.)
