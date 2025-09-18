# Model_outputs
# TARSS model outputs

This repository contains code and **model outputs** (tables/figures) for our analyses (2014–2022), including MDR.

## Contents
- `code/`: R scripts to reproduce summaries/plots.
- `results/tables`: Excel exports of model summaries and comparators.

## Reproducibility
- R ≥ 4.2, packages: `tidyverse`, `broom`, `janitor`, `readxl`, `writexl`, `ggplot2`, `car`, `emmeans`.
- Run scripts in order: `01_data_prep.R` → `02_models_mdr.R` → `03_export_tables_figs.R`.

## Licensing
- Code: MIT. Outputs: CC BY 4.0.
