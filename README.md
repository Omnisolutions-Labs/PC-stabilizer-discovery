# PC-stabilizer-discovery

Chemoinformatics-guided discovery of food-grade anionic stabilizers for phycocyanin under acidic conditions.

This repository accompanies the manuscript:

> Chuang, K.; Luo, L.; Law, L. *Chemoinformatics-guided discovery of food-grade anionic stabilizers for phycocyanin under acidic conditions: integrating QSPR, a chemistry-prior heuristic, and virtual screening of GRAS additives.* (submitted, 2026).

A 48-compound training dataset, a LightGBM + RDKit QSPR pipeline, a chemistry-prior heuristic baseline, virtual screening of 30 GRAS food additives, and 24-bottle wet-lab validation data are all provided, along with the figures and code used to produce them.

[![License: MIT](https://img.shields.io/badge/Code%20License-MIT-blue.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/Data%20License-CC%20BY%204.0-lightgrey.svg)](LICENSE-DATA)
[![DOI](https://img.shields.io/badge/Zenodo%20DOI-pending-orange.svg)](https://zenodo.org/)

## Highlights

- **First QSPR-ready open dataset** for phycocyanin acid-stabilizer discovery (48 unique compounds, 3 effective hits).
- Descriptor-based **LightGBM LOO-AUC = 0.73** (marginally above a single-feature baseline at 0.67); a chemistry-prior heuristic encoding anion-type priors reaches **AUC = 0.95** — illustrating that, in this genuinely small-data regime, explicit chemical knowledge outperforms descriptor-based ML.
- **Virtual screening of 30 GRAS food additives** ranked the pyrophosphate family at the top.
- **Wet-lab validation (n = 3 × 8 groups)** confirmed sodium hexametaphosphate (SHMP, 78.1% color retention), tetrasodium pyrophosphate (TSPP, 54.1%), and sodium phytate (IP₆, 52.7%) as new effective stabilizers; an IP₆ + STPP ternary combination reached 83.8% retention.
- **Mechanistic evidence**: per-bottle ζ-potential correlates with color retention (Pearson r = −0.82, p = 1 × 10⁻⁶, n = 24).

## Repository structure

```
PC-stabilizer-discovery/
├── README.md                            ← this file
├── LICENSE                              ← MIT (code)
├── LICENSE-DATA                         ← CC BY 4.0 (data)
├── CITATION.cff                         ← citation metadata
├── requirements.txt                     ← Python dependencies
├── .gitignore
│
├── data/
│   ├── training/
│   │   ├── curated_compounds.csv        ← 48-compound dataset with SMILES + labels
│   │   ├── features_matrix.csv          ← 21-dimensional feature matrix
│   │   ├── training_predictions.csv     ← QSPR + chemistry-prior probabilities (LOO-CV)
│   │   └── model_metrics.json           ← AUC, accuracy, n
│   ├── virtual_screen/
│   │   └── virtual_screening_results.csv  ← 30 GRAS candidates ranked by QSPR
│   └── validation/
│       ├── UV_long.csv                  ← per-bottle UV-Vis raw data (24 bottles × 2 timepoints)
│       ├── Zeta_final.csv               ← ζ-potential measurements (final, with rerun handling)
│       ├── summary.csv                  ← group-level mean ± SD (Table 4 in the paper)
│       ├── predicted_vs_observed.csv    ← Table 4 prediction vs experiment comparison
│       └── ttest_vs_G0.csv              ← Welch's t-test results vs negative control
│
├── code/
│   ├── qspr_pipeline.py                 ← end-to-end QSPR training, SHAP, virtual screening
│   └── analyze_real.py                  ← validation-data analysis + Figure 6 generation
│
├── figures/
│   ├── Fig1_dataset_overview.png        ← 3-round iterative dataset
│   ├── Fig2_SHAP_summary.png            ← feature importance
│   ├── Fig3_ROC_comparison.png          ← QSPR vs chemistry-prior vs baseline
│   ├── Fig4_virtual_screen.png          ← GRAS virtual screening top-12
│   ├── Fig5_QSPR_vs_chemprior.png       ← cross-method agreement scatter
│   └── Figure6_validation_panel.png     ← wet-lab validation (4 sub-panels)
│
├── sop/
│   └── validation_experiment_SOP.md     ← experimental protocol (preregistration evidence)
│
└── docs/
    └── methods_supplement.md            ← extended methods, notes, parameter rationale
```

## Quick start

### 1. Reproduce the QSPR analysis

```bash
git clone https://github.com/Omnisolutions-Labs/PC-stabilizer-discovery.git
cd PC-stabilizer-discovery
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
python code/qspr_pipeline.py
```

Outputs:
- LOO-CV AUC for the three methods (QSPR LightGBM, chemistry-prior heuristic, single-feature baseline)
- SHAP feature importance figure
- Virtual-screening rank list for the 30-compound GRAS library

Random seed is fixed (`SEED = 20260519`); all results are bit-level reproducible.

### 2. Reproduce the validation analysis and Figure 6

```bash
python code/analyze_real.py
```

This reads `data/validation/UV_long.csv` and `data/validation/Zeta_final.csv`, computes per-bottle CR%, group-level statistics, the Pearson correlation between ζ-potential and color retention, and regenerates the 4-panel Figure 6.

## Data dictionary

### `data/training/curated_compounds.csv`

| Column | Description |
|---|---|
| `compound_name_en` | English compound name |
| `compound_name_cn` | Chinese compound name |
| `category` | Chemical class (salt, polysaccharide, polyphosphate, ...) |
| `SMILES` | SMILES string (empty for polymers without canonical structure) |
| `mw` | Molecular weight (g/mol; estimated for polymers) |
| `is_polymer` | 0/1 flag |
| `anion_type` | Anionic functional group class |
| `functional_group_density` | Anionic groups per repeat unit (or per molecule) |
| `best_observed_label` | 1 = effective; 0 = failed; −1 = precipitate |
| `source_round` | Screening round identifier |
| `notes` | Laboratory observations |

### `data/validation/UV_long.csv`

Per-bottle UV-Vis data for the 24-bottle validation experiment (8 groups × n = 3) at t = 0 and t = 7 d, including A₆₂₀, A₃₅₀, A₇₀₀, blank values, computed color retention, raw file references, and operator/timestamp metadata.

### `data/validation/Zeta_final.csv`

Per-bottle ζ-potential measurements (three within-instrument runs each), polydispersity index, derived intensity, and notes about excluded measurements (B14, B21 re-acquired per SOP).

## Citation

If you use this dataset, code, or framework in your work, please cite:

```bibtex
@article{Chuang2026PCStabilizer,
  author  = {Chuang, Kevin and Luo, Lorry and Law, Luke},
  title   = {Chemoinformatics-guided discovery of food-grade anionic stabilizers for phycocyanin under acidic conditions: integrating QSPR, a chemistry-prior heuristic, and virtual screening of GRAS additives},
  journal = {ACS Food Science \& Technology},
  year    = {2026},
  note    = {Manuscript submitted},
}
```

Once a Zenodo DOI is minted, please cite both the article and the dataset DOI.

## Licenses

- **Code** (everything under `code/`): MIT License — see [LICENSE](LICENSE).
- **Data and figures** (everything under `data/`, `figures/`, `sop/`, `docs/`): Creative Commons Attribution 4.0 International — see [LICENSE-DATA](LICENSE-DATA).

Both licenses permit unrestricted academic and commercial reuse with attribution.

## Funding and conflict of interest

This research was funded by **Omnisolutions Laboratory Holdings Limited** (Hong Kong SAR, China). All authors are employees of the funder.

## Contact

For questions about the dataset or code, please contact **Luke Law** at <luke.law@omnisolutionslabs.com>.

## Acknowledgments

We thank lab technicians LX and ZHY for diligent execution of the validation experiments.
