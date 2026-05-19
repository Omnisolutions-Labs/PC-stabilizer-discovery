# Methods supplement

This document expands on the Methods section of the manuscript with implementation details, parameter rationales, and reproducibility notes.

## A. Feature engineering

The 21-dimensional feature vector combines 10 RDKit descriptors and 11 hand-curated domain-expert features.

### A.1 RDKit descriptors (10)

| Feature | RDKit call | Notes |
|---|---|---|
| MolWt | `Descriptors.MolWt(mol)` | g/mol |
| LogP | `Crippen.MolLogP(mol)` | dimensionless |
| TPSA | `Descriptors.TPSA(mol)` | Å² |
| HBD | `Lipinski.NumHDonors(mol)` | count |
| HBA | `Lipinski.NumHAcceptors(mol)` | count |
| FormalCharge | sum of `GetFormalCharge()` over atoms | integer |
| NumRings | `Lipinski.RingCount(mol)` | count |
| NumRotBonds | `Lipinski.NumRotatableBonds(mol)` | count |
| FracSP3 | `rdMolDescriptors.CalcFractionCSP3(mol)` | 0–1 |
| NumHeavyAtoms | `mol.GetNumHeavyAtoms()` | count |

Polymer compounds without a canonical SMILES are assigned default RDKit values:
```python
MolWt=mw, LogP=-2.0, TPSA=200, HBD=10, HBA=15,
FormalCharge=-int(eff_neg_charge_per_unit*5),
NumRings=1, NumRotBonds=20, FracSP3=0.9, NumHeavyAtoms=50
```

### A.2 Domain-expert features (11)

| Feature | Definition |
|---|---|
| is_polymer | 0/1 |
| log10_mw | log₁₀ of molecular weight |
| is_phosphate | 1 if anion class ∈ {phosphate, polyphosphate, metaphosphate} |
| is_polyphosphate | 1 if anion class = polyphosphate |
| is_sulfate / is_sulfonate / is_carboxylate | one-hot encoding |
| is_anionic | 1 if any anionic group present |
| fg_density | functional groups per repeat unit |
| eff_neg_charge_per_unit | fg_density × pKa-weighted ionization fraction at pH 3 |
| is_metal_polyvalent | 1 if compound contains Ca²⁺, Mg²⁺, Zn²⁺, Fe³⁺, etc. |

The effective negative charge weighting at pH 3 uses:
- Polyphosphate, sulfate, sulfonate: 1.0 (fully ionized)
- Metaphosphate: 0.7
- Carboxylate: 0.4 (partial ionization; pKa ≈ 3–4)
- Hydroxyl, amide, none: 0.0

## B. Model hyperparameters and training

LightGBM classifier:
```python
LGBMClassifier(
    n_estimators=200,
    learning_rate=0.04,
    num_leaves=15,
    min_child_samples=2,
    max_depth=4,
    reg_lambda=1.0,
    class_weight={0: 1.0, 1: pos_weight},  # pos_weight = n_neg / n_pos
    random_state=20260519
)
```

`pos_weight` is computed dynamically as `(n_total − n_positive) / n_positive` ≈ 11.25 for the 49-compound dataset. Leave-one-out cross-validation is used due to the small dataset size; 49 LightGBM models are trained, each on 48 examples.

## C. SHAP attribution

`shap.TreeExplainer` is applied to the final model (trained on all 49 examples) to obtain Shapley values for each feature on each compound. The summary plot reflects the absolute Shapley contributions across the full dataset.

## D. Chemistry-prior heuristic

The chemistry-prior heuristic is a deterministic rule set assigning P(effective) directly from compound name and anion class (see `code/qspr_pipeline.py`, function `chem_prior` / `llm_like_zero_shot`). It encodes published pKa values and reported activity of major anion classes at pH 3.

Probabilities are calibrated against the 49-compound training set; refinement against an external benchmark dataset is left to future work.

## E. Virtual screening library

Thirty GRAS food additives were curated by combining EFSA E-number lists with the Chinese GB 2760 additive catalog, restricted to anionic polymers and phosphates. Polymers without canonical SMILES use the polymer default descriptor set described in §A.1.

## F. Statistical analysis details

- One-sided Welch's t-test (`scipy.stats.ttest_ind(..., equal_var=False, alternative='greater')`) for treatment vs G0 contrasts.
- Pearson r and two-tailed p-value via `scipy.stats.pearsonr` for ζ vs CR₆₂₀.
- Holm-Bonferroni correction across the seven vs-G0 contrasts at family-wise α = 0.05.

## G. Reproducibility

Random seed: `SEED = 20260519`. All NumPy random generators, LightGBM `random_state`, and Monte-Carlo sampling routines use this seed. Re-running `python code/qspr_pipeline.py` and `python code/analyze_real.py` reproduces all reported AUC, accuracy, virtual-screening ranks, and Figure 6 panels at bit-level (verified on Python 3.13, numpy 1.26, pandas 2.0, scipy 1.15.3, scikit-learn 1.5, lightgbm 4.0, shap 0.45).

## H. Limitations

1. Class imbalance (4 positives / 45 negatives) limits LOO-CV power; SHAP attributions are interpreted qualitatively.
2. Polymer compounds with default descriptors may underestimate within-class differences.
3. Sodium phytate (IP₆) was a known descriptor-model false negative; the model's applicability domain is therefore limited to linear polyphosphates and acyclic polyanions.
4. The 7-day accelerated assay at 46 °C is a surrogate for the industrially relevant 30-day shelf life at 25 °C; real-product testing is outside the present scope.
