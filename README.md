# Water Quality Analysis for 15 Wells

This repository contains an R workflow for analyzing water-quality data from 15 wells. The script reads an Excel workbook, cleans the dataset, calculates descriptive statistics, checks compliance against water-quality standards, performs inferential statistics, and generates charts and export files.

## Project purpose

The analysis evaluates several water-quality indicators, including:

- pH
- Total dissolved solids (TDS)
- Electrical conductivity (EC)
- Turbidity
- Total coliform count (TCC)
- Nitrate
- E. coli
- Distance from the septic tank

The workflow is implemented in [kamal.R](kamal.R) and produces outputs in the [output](output) folder.

## Files

- [kamal.R](kamal.R): Main R script that performs the analysis and exports results.
- [kamal.xlsx](kamal.xlsx): Source dataset in Excel format.
- [output](output): Folder containing generated Excel reports and plots.

## Dependencies

The script uses the following R packages:

- readxl
- dplyr
- tidyr
- ggplot2
- writexl

Install them in R with:

```r
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "writexl"))
```

## How to run

1. Open R or RStudio in the repository root.
2. Run the analysis script:

```r
source("kamal.R")
```

3. The script will create the output folder if needed and save all generated files there.

## What the script produces

The analysis generates:

- Descriptive statistics in Excel
- Compliance checks against WHO/NSDWQ-style criteria
- A comparison table relating distance to contamination indicators
- Correlation analyses (Pearson and Spearman)
- Normality testing (Shapiro-Wilk)
- Simple linear regression results for distance versus contaminants
- Scatter plots for key relationships
- Group comparisons between close and far wells
- A minimum safe distance assessment

## Output folder

The generated files are written to [output](output), including:

- Excel files with descriptive and inferential statistics
- PNG plots for visual assessment

## Notes

- The script assumes that the input workbook contains the expected water-quality columns and that the first empty column is removed during preprocessing.
- Results are intended for exploratory analysis and reporting, and should be interpreted with the sample size and context of the study in mind.
