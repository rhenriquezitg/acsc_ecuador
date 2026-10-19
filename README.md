# Hospitalizations for ambulatory care sensitive conditions in Ecuador, 2000–2024

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23366479.svg)](https://doi.org/10.5281/zenodo.23366479)

Data, code and supplementary material for:

> Martínez-Abarca AL, Henríquez-Trujillo MC, Bonilla E, Henríquez-Trujillo AR. Hospitalizations for ambulatory care sensitive conditions in Ecuador, 2000–2024. Submitted to Revista de Saúde Pública.

The study describes hospitalizations for ambulatory care sensitive conditions (ACSC) in Ecuador from 2000 to 2024. It uses all public and private hospital discharges recorded by the Instituto Nacional de Estadística y Censos (INEC) and examines:

- national trends (joinpoint regression);
- changes around the 2008 health reform and the onset of the COVID-19 pandemic (segmented regression);
- geographic variation between the 24 provinces;
- the association between hospital bed supply and hospitalizations by province, 2000–2024.

## Contents

| Path | Content |
|---|---|
| `acsc_ecuador_analysis.R` | Single script that reproduces every number, table and figure of the article and its supplementary material |
| `supplementary/Supplementary_material.docx` (and `.pdf`) | Tables S1–S3 and Figures S1–S8 |
| `supplementary/Table_S4.xlsx` | Joinpoint regression of all series |
| `supplementary/Table_S5.xlsx` | Segmented regression by ACSC group |
| `supplementary/Table_S6.xlsx` | Sensitivity analyses |
| `supplementary/Table_S7.xlsx` | Annual age-standardized rates by province, 2000–2024, with 95% confidence intervals |
| `supplementary/Table_S8.xlsx` | Annual indicators of geographic variation |
| `supplementary/Table_S9.xlsx` | Hospital bed supply by province and its association with hospitalizations, 2000–2024 |
| `data/` | Aggregated input data (see below) |
| `results/` | All results produced by the script, as CSV files |
| `figures/` | Figures 1–2 of the article and Figures S1–S8 (TIFF 300 dpi and PNG) |

## Data

`data/` holds aggregated counts, not record-level data:

| File | Description |
|---|---|
| `acsc_annual.csv.gz` | Hospitalizations by place basis (`geo`: `res` residence, `occ` occurrence), year, province, sex, single year of age and group |
| `acsc_monthly.csv.gz` | National hospitalizations by year, month of discharge, age and group |
| `acsc_icd.csv.gz` | ACSC hospitalizations by year, province of residence, sex, age, group and ICD-10 code of the main diagnosis |
| `record_log.csv` | Records with missing or invalid values, by year |
| `population.csv.gz` | INEC population projections (revision 2024) by year, province, sex and single year of age |
| `standard_population.csv` | Standard populations: Ecuador 2022 projected population (main analysis) and WHO world standard (sensitivity) |
| `icd10_labels.csv` | ICD-10 labels for the codes used |
| `beds_province.csv` | Available beds and beds of normal allocation by province of the establishment, 2000–2024: INEC establishment files for 2000–2003 and series 1.5.2 and 1.5.3 for 2004–2024 |
| `beds_unassigned_2000_2003.csv` | Beds of establishments with an unassigned province code in the 2001–2003 files (left out) |
| `beds_national.csv` | National population, beds, discharges and bed occupancy, 1997–2024 (INEC series 1.4.2) |
| `provinces_geo.rds` | Provincial boundaries (ADM1 2024) and neighbour lists (queen contiguity and 4 nearest neighbours) |

Group codes: 1–20 are the ACSC groups (defined in section 1 of the script); 99 is deliveries (O80–O84); 0 is all other causes. Sex: 1 men, 2 women. Provinces use the INEC codes 01–24 in current boundaries; 88 and 90 are abroad and non-delimited zones. Age is in completed years, with 100 meaning 100 or older.

Bed counts come from INEC's hospital beds registry: the establishment-level files for 2000–2003 (`bdd_camas_hosp_<year>_spss.zip`, rebuilt with `beds_2000_2003()` in section 3) and the published historical series for 2004–2024 (*Tabulados series ECEH 2024*). The source microdata are the INEC hospital discharge files (*Egresos Hospitalarios*, later *Camas y Egresos Hospitalarios*), 2000–2024. They are publicly available at [ecuadorencifras.gob.ec](https://www.ecuadorencifras.gob.ec). The authors harmonized the files across years and aggregated them with an R package of their own, `inectab`. Section 3 of the script rebuilds `data/` from the harmonized microdata when `REBUILD_FROM_MICRODATA <- TRUE`. The rest of the analysis needs only the files in `data/`.

## Methods in brief

- **ACSC list:** 20 groups of the Brazilian list (Alfradique et al. 2009), classified on the main diagnosis.
- **Rates:** directly age-standardized, with 18 five-year age bands and Ecuador's 2022 projected population as the standard. Confidence intervals by the gamma method (Fay and Feuer).
- **Joinpoint regression:** weighted log-linear models. The number of joinpoints is chosen by permutation tests (Kim et al. 2000).
- **Segmented regression:** quasi-Poisson models with Fourier seasonal terms and Newey–West standard errors.
  - 2008 reform: January 2000–December 2019, interruption in January 2008.
  - COVID-19: January 2012–December 2024, interruption in March 2020.
- **Hospital beds:** annual provincial counts, 2000–2024, modelled with quasi-Poisson regression with province and year fixed effects, log beds per 10,000 as exposure and expected counts as offset; standard errors clustered by province. Bed counts are smoothed with a repeated five-year running median to remove registration spikes. For 2000–2003, Santo Domingo de los Tsáchilas and Santa Elena (created in 2007) receive their 2004 share of the beds of Pichincha and Guayas.
- **Geographic variation:** coefficient of variation, extremal quotient, systematic component of variation, and global and local Moran's I.

## How to reproduce

R 4.2 or later. Install the packages once, then run the script from the repository root:

```r
install.packages(c("dplyr", "tidyr", "readr", "ggplot2", "patchwork", "sf", "spdep",
                   "sandwich", "ragg", "scales"))
```

```sh
Rscript acsc_ecuador_analysis.R
```

The script takes a few minutes and writes `results/` and `figures/`. The package versions used are listed in `results/session_info.txt`.

## License

- Code: MIT License (`LICENSE`).
- Aggregated data, results, figures and supplementary material: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
- Source data: INEC, Ecuador.

## Citation

Please cite the article and this repository: https://doi.org/10.5281/zenodo.23366479 (see `CITATION.cff`).

## Contact

Aquiles Rodrigo Henríquez-Trujillo, Institute of Tropical Medicine, Antwerp. henriquezmd@gmail.com. ORCID [0000-0002-3094-4438](https://orcid.org/0000-0002-3094-4438).
