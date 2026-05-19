# Validation Experiment Standard Operating Procedure (SOP)

> Companion to: *Chemoinformatics-guided discovery of food-grade anionic stabilizers for phycocyanin under acidic conditions* (Chuang et al., 2026).

This SOP corresponds to the 24-bottle wet-lab validation experiment reported in Section 3.6 of the manuscript. All quantitative data from this protocol are released in `data/validation/`.

## 1. Scope and acceptance criteria

Validate the QSPR- and chemistry-prior-predicted top candidates (SHMP, TSPP, sodium phytate) at pH 3 and 46 °C over 7 days.

Acceptance criteria (any one yields a primary result):
- SHMP (V1-a or V1-b) CR₆₂₀ ≥ 60%
- TSPP (V2-a or V2-b) CR₆₂₀ ≥ 50%
- Sodium phytate (V3) CR₆₂₀ ≥ 50%
- ζ-potential of any treatment group at pH 3 below −10 mV

## 2. Experimental matrix (24 bottles, n = 3)

| Group | Treatment | Ratio (stab:PC) | n | Bottle IDs |
|---|---|---|---|---|
| G0     | PC + HCl (negative)        | —      | 3 | B01–B03 |
| G_pos  | STPP (positive control)    | 2:1    | 3 | B04–B06 |
| V1-a   | SHMP (top-1 QSPR)          | 1:1    | 3 | B07–B09 |
| V1-b   | SHMP                       | 2:1    | 3 | B10–B12 |
| V2-a   | TSPP                       | 1:1    | 3 | B13–B15 |
| V2-b   | TSPP                       | 2:1    | 3 | B16–B18 |
| V3     | Sodium phytate (IP₆)       | 1:1    | 3 | B19–B21 |
| V3-pos | IP₆ + STPP (ternary)       | 1:1:1  | 3 | B22–B24 |

## 3. Critical control points (CCP)

| CCP | Parameter | Locked value | Failure mode if violated |
|---|---|---|---|
| CCP-1 | PC final concentration | 0.375 % w/w | Baseline drift; CR not comparable |
| CCP-2 | Target pH | 3.00 ± 0.05 | ζ-titration inflection shifts ± 5 mV |
| CCP-3 | Incubation temp × time | 46.0 ± 0.5 °C × 168 h | Kinetics drift |
| CCP-4 | Sample / headspace | 10.0 mL / 14 mL | Oxidation accelerates fading |
| CCP-5 | PC stock thaw/refreeze | Single use, same day | Partial denaturation |
| CCP-6 | HCl titration rate | 1.00 M, 1 µL/s, 200 rpm stir | Local pH < 2 precipitation |

## 4. Materials

Food-grade phycocyanin E18+ (Yunnan Lü-A Bio, lot PC-E18-260509-B), STPP (Aladdin, lot STPP-260428-03), SHMP (Macklin, lot SHMP-260502-11), TSPP (Macklin, lot TSPP-260430-07), sodium phytate (Sigma, lot IP6-260421-02), 1.00 M HCl (titrated), 24 mL amber headspace vials, 0.22 µm PVDF filters, NIST-traceable pH calibration buffers.

## 5. Workflow (10-day timeline)

- **D-2** Materials ready; vials cleaned, dried, and labelled (B01–B24).
- **D-1** Stocks prepared: 0.50% w/w PC, 75 mg/mL stabilizer master solutions, 1 mM NaCl pH 3 ζ-dilution buffer; HOBO logger placed in incubator 24 h preheat.
- **D0**  Sample preparation, t = 0 UV-Vis scan, N₂ purge, seal, incubate.
- **D1–D6** Twice-daily temperature check (09:00 / 17:00); HOBO logger continuous.
- **D7** End-point: ambient equilibration 30 min → UV-Vis full scan → D65 light-box photo → ζ-potential (1:10 dilution into 1 mM NaCl pH 3; 3 within-instrument repeats; latex standard every 5 samples).
- **D8** Data entry, blank correction, CR₆₂₀ calculation.
- **D9** Statistical analysis (Welch's t-tests vs G0; Pearson ζ-CR correlation).

## 6. Sample preparation (per bottle, 10.0 mL final volume)

| Group | PC master (0.50%) | Stabilizer master (10×) | 1 M HCl | ddH₂O to |
|---|---|---|---|---|
| G0     | 7.50 mL | —                                  | ≈ 0.10 mL | 10.00 mL |
| G_pos  | 7.50 mL | 1.00 mL STPP (75 mg/mL)            | ≈ 0.13 mL | 10.00 mL |
| V1-a   | 7.50 mL | 0.50 mL SHMP                       | ≈ 0.13 mL | 10.00 mL |
| V1-b   | 7.50 mL | 1.00 mL SHMP                       | ≈ 0.13 mL | 10.00 mL |
| V2-a   | 7.50 mL | 0.50 mL TSPP                       | ≈ 0.13 mL | 10.00 mL |
| V2-b   | 7.50 mL | 1.00 mL TSPP                       | ≈ 0.13 mL | 10.00 mL |
| V3     | 7.50 mL | 1.00 mL IP₆ (37.5 mg/mL)           | ≈ 0.15 mL | 10.00 mL |
| V3-pos | 7.50 mL | 0.5 mL IP₆ + 0.5 mL STPP           | ≈ 0.15 mL | 10.00 mL |

## 7. Measurement protocols

**UV-Vis** — 1 cm quartz cell; spectrum 350–750 nm, step 1 nm, medium speed; blank = matching pH and stabilizer concentration without PC. CR₆₂₀ = [A₆₂₀(7d) − A_blank] / [A₆₂₀(0) − A_blank] × 100.

**ζ-Potential** — Malvern Zetasizer Nano ZS with DTS1070 folded-capillary cells; samples diluted 1:10 (or 1:20 if PDI > 0.5) in 1 mM NaCl, pH 3.00; 25 °C; 60 s equilibration; three runs per bottle. Latex standard target −55 ± 5 mV verified every five samples.

**Photo documentation** — D65 light box, ISO 100, f/8, 1/125 s, 50 mm; gray-card white balance; Pantone 285/286/2935 blue color cards as reference.

## 8. Statistical analysis

- One-sided Welch's t-test vs G0 for each treatment (H₁: treatment > G0)
- Pearson r between per-bottle ζ-potential and CR₆₂₀
- Holm-Bonferroni adjustment across the seven vs-G0 contrasts at family-wise α = 0.05
- Synergy test: one-sided Welch's t-test of V3-pos vs V3 alone

## 9. Pre-registration

This SOP was finalized and recorded internally on 2026-05-11 prior to the start of bottle preparation (2026-05-12), serving as a preregistration document for the validation campaign.
