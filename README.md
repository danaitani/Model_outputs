# Model_outputs
# TARSS model outputs

This repo contains the analysis for pathogen–antibiotic resistance and MDR using multivariable logistic regression in R.

/Outcomes/
Resistance: EUCAST interpretation, modelled as R vs S/I.
MDR: per Magiorakos et al., 2012 — non-susceptibility (R/I) to ≥1 agent in ≥3 classes.

/Core covariates & references/
Year (continuous, per 1-year increase)
Specimen type (ref = Blood)
Hospital (ref = HCN)
Collapsed ward (Ward_ED_ICU: ref = ICU; levels: ICU / ED / Other)
Gender (ref = Male) — used in sensitivity models
Age (7-level) (ref = 30–<52 y; bands: 0–28 d, 29–365 d, 1–5 y, 6–<30 y, 30–<52 y, 52–<67 y, ≥67 y) — used in sensitivity models
Ward categories were collapsed (ICU, ED, Other) based on counts; see the collapsed-check table in the results folder.

/Model set (additive → interaction → +demographics)/
Model	Outcome	Formula (schematic)	Description / notes
M1	Resistance	Resistance ~ Year + Specimen_type + Hospital + Ward_ED_ICU	Additive main effects (collapsed ward: ICU / ED / Other).
M1′	MDR	MDR ~ Year + Specimen_type + Hospital + Ward_ED_ICU	Additive MDR analogue.
M2	Resistance	Resistance ~ Year + Specimen_type + Hospital * Ward_ED_ICU	Primary interaction model (hospital×ward).
M2′	MDR	MDR ~ Year + Specimen_type + Hospital * Ward_ED_ICU	MDR interaction analogue.
M3	Resistance	M2 + Gender	Sensitivity: adds Gender (ref: Male).
M3′	MDR	M2′ + Gender	MDR analogue.
M4	Resistance	M2 + Age(7-level)	Sensitivity: adds Age (ref: 30–<52 y); exploratory due to <1 y missingness.
M4′	MDR	M2′ + Age(7-level)	MDR analogue.

/Model selection & diagnostics/
Nested comparisons (e.g., M1 vs M2; M2 vs M3/M4): Likelihood-ratio tests (LRT) on the same rows; ΔAIC > 2 considered meaningful.
Discrimination: AUC (excellent > 0.90; acceptable 0.80–0.90; fair 0.70–0.80; poor < 0.70).
Collinearity: VIF (car), flag if VIF > 5.
CIs: Wald 95% CIs (robust when profile CIs fail under separation).

/Sensitivity analysis/
We evaluated Gender and Age (7-level) on top of the interaction model (M2/M2′). Only covariates that improved fit (LRT p < 0.05 and/or ΔAIC > 2) were retained. 

Reproducibility
R (≥ 4.x); key packages: dplyr, tidyr, broom, car, pROC, writexl, forcats, AMR.
ORs and CIs exported via broom::tidy() (Wald CIs).
Outputs (ORs, global p-values, AIC/AUC, VIF, comparators) are written to /results/... as .xlsx.
