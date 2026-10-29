# dataMineR

`dataMineR` provides clinical data-mining and comparison tools for R. The
package supports grouped Table 1 summaries, regression screening, continuous
and categorical group comparisons, correlation and agreement analyses,
survival and mediation analyses, dimensionality reduction, power analysis,
and bootstrap predictor-stability screening.

## Installation

```r
# install.packages("remotes")
remotes::install_github("brendangongol/dataMineR")

# Install optional packages used by the complete demonstration workflow:
# remotes::install_github("brendangongol/dataMineR", dependencies = TRUE)
```

## Getting started

```r
library(dataMineR)

load(system.file("extdata", "dfMigrant.rda", package = "dataMineR"))
InitialTablebg(dfMigrant, GroupBy = "Migrant", p_values = TRUE)
```

The bundled workflow demonstrates Table 1 generation, normality checks,
group comparisons, regression screening, survival and mediation analyses,
dimensionality reduction, and power analysis with the synthetic `dfMigrant`
dataset:

```r
demo("dataMineRWorkflow", package = "dataMineR")
```

or open [`demo/dataMineRWorkflow.R`](demo/dataMineRWorkflow.R) directly.
The demo writes to `dataMineRWorkflow-results/` under the current directory.
Set `DATAMINER_OUTPUT_DIR` and `DATAMINER_WORKFLOW_NAME` to use a different
root or keep multiple workflow outputs separate. Set `DATAMINER_SEED` to
control the workflow's random seed (the default is `2026`).

## Example datasets

The package ships with several synthetic datasets used in the demo and
function examples: `dfMigrant`, `copd_data`, `dream_data`, `peru_data`, and
`cardiac_data`. It also includes a lightly reshaped copy of
`medicaldata::smartpill` (`MData`) and an NHANES adult subset
(`nhanes_like`). See `data-raw/create_example_data.R` for the generation
code.

## Development

Functions are organized by topic in `R/`. Documentation is maintained in
roxygen2 comments above each function. After changing those comments,
regenerate `NAMESPACE` and `man/` with:

```r
roxygen2::roxygenise()
devtools::test()
devtools::check()
```
