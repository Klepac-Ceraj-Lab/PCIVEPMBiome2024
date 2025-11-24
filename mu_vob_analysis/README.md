# MU_VOB - Reproducible Stats

This folder contains code to reproduce the analyses for the paper.

## Contents

- `replicate.py` - Python analysis script
- `replicate.R` - R analysis script (optional alternative)
- `requirements.txt` - Python package versions

## Input data

Input data files should be sourced/download from the
centralized data repository in Dada Dryad (DOI 10.5061/dryad.15dv41p9k).
The following data files are expected in the working directory:

- `MU_VOB.csv` - dataset
- `MU_VOB_codebook.csv` - machine-generated codebook
- `column_manifest.json` - auto-detected column names used by scripts

## Quickstart (Python)

```bash
python3 -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt

# run on the included CSV
python replicate.py MU_VOB.csv

# or on a different CSV with the same columns
python replicate.py path/to/your.csv
```

The script writes CSV tables to `results/`:
- `correlations_r.csv`, `correlations_p.csv`
- `ANCOVA_typeII.csv`
- `AgeOnly_ANOVA.csv`, `ResidualCorr.txt`
- `ChiSquare_counts.csv`, `ChiSquare_results.csv`
- `IC_TypeIII_ANOVA.csv`
- `IC_simple_slopes_raw.csv`, `IC_simple_slopes_standardized.csv`

## Notes

- Correlations use series-mean imputation (matching SPSS “replace missing with series mean”).  
- The inhibitory-control model uses Type III sums of squares and the predictors:
  `child_sex`, `epds`, `VOB_For_Age_Residual`, `age_pci_days_NoMissingValues`, `MU`,
  plus `age×MU`, `VOBres×age`, `VOBres×MU`, and `VOBres×age×MU`.
- Simple slopes of MU are evaluated at ±1 SD of age and VOB-for-age with HC3 robust SEs.
- If your column names differ slightly (e.g., `VOB_raw` instead of `VOB`), the script
  reads `column_manifest.json` to find the correct columns. Edit it if needed.

Generated on 2025-11-14.
