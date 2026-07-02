# Computational Benchmarking of the INLA-SPDE Approach for Spatial Point Processes: Parameter Recovery and Spatial Confounding

**Robyn Irawan**  
Center for Mathematics and Society, Faculty of Science,  
Parahyangan Catholic University (UNPAR), Bandung, West Java, Indonesia  
robynirawan.tjia@unpar.ac.id

---

## Overview

This repository contains the complete R pipeline for the Monte Carlo simulation study reported in:

> Irawan, R. (2026). Computational Benchmarking of the INLA-SPDE Approach for Spatial Point Processes: Parameter Recovery and Spatial Confounding. *Journal of Statistical Computation and Simulation* (under revision, Manuscript ID: GSCS-2026-0551).

The study benchmarks the inferential performance of the Log-Gaussian Cox Process (LGCP) estimated via the INLA-SPDE framework, with focus on parameter recovery, spatial confounding, and spatial range identifiability on bounded domains. A 500-replicate Monte Carlo simulation is conducted on a synthetic spatial domain with known ground-truth parameters, enabling controlled evaluation of empirical credible interval coverage, bias, and RMSE for fixed-effect covariates and Matern spatial hyperparameters.

---

## Repository Contents

```
lgcp-inla-spde-benchmarking/
├── README.md
├── LGCP_INLA_SPDE_Benchmarking.Rmd        # Complete R pipeline
└── outputs/
    └── tables/
        ├── 00_covariate_standardisation_stats.csv
        ├── 00_mesh_diagnostics.csv
        ├── 00_reference_encounter_counts.csv
        ├── 01_reference_fixed_effects.csv
        ├── 02_reference_hyperparameters.csv
        ├── 03_reference_recovery.csv
        ├── 04_reference_waic_logcpo.csv
        ├── 05_reference_delta_waic.csv
        ├── 06_pc_prior_sensitivity.csv
        ├── 07_simulation_summary.csv
        ├── 08_simulation_summary_with_mcse.csv
        └── 13_appendix_A_reference_coordinates.csv
```

The `outputs/tables/` folder contains pre-computed results from the full 500-replicate run, allowing inspection of all reported numerical results without re-running the full pipeline. Running the `.Rmd` from scratch will regenerate all tables and figures.

---

## Software Requirements

| Software  | Version used   |
|-----------|----------------|
| R         | 4.4.0          |
| INLA      | 24.06.27       |
| sf        | 1.1-0          |
| spatstat  | 3.5-1          |
| ggplot2   | CRAN current   |
| patchwork | CRAN current   |
| viridis   | CRAN current   |
| Matrix    | CRAN current   |
| fields    | CRAN current   |

### Installing R-INLA

INLA is not on CRAN and must be installed from the INLA repository:

```r
install.packages("INLA",
  repos = c(getOption("repos"),
            INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE, type = "binary")
```

All other packages can be installed via:

```r
install.packages(c("sf", "spatstat", "spatstat.geom", "spatstat.random",
                   "ggplot2", "patchwork", "viridis", "Matrix", "fields"))
```

The `.Rmd` file contains an `installation` chunk (set to `eval=FALSE`) with the full installation code. Run that chunk once before the first use.

---

## How to Run

1. Clone or download this repository.
2. Open `LGCP_INLA_SPDE_Benchmarking.Rmd` in RStudio.
3. Set the working directory to the project root folder via Session → Set Working Directory → To Source File Location.
4. Run all chunks sequentially, or click **Knit** to produce a full PDF report.

Output folders are created automatically on first run:

```
outputs/tables/       <- CSV result tables
outputs/figures/      <- PDF figures
outputs/checkpoints/  <- simulation state saved every 25 replicates
```

### Expected Runtime

| Component                        | Approximate time |
|----------------------------------|-----------------|
| Reference run and diagnostics    | 5–10 minutes    |
| 500-replicate Monte Carlo study  | 5–7 hours       |
| Figures and final summary        | 2–5 minutes     |

Tested on: Windows 11 x64, AMD Ryzen 7 8840HS, 16 GB RAM.

Checkpoints are saved every 25 replicates so the simulation can be safely interrupted and resumed.

---

## Simulation Design

**Spatial domain:** Synthetic 10 × 10 unit square with three ecologically distinct habitat types (Water, Rock/Cave, Grass), inspired by the Kanto region topology (Game Freak, 1996).

**Data-generating process:** Single-covariate LGCP,

```
log λ_k(s) = β0,k + β1,k · X̃_k(s) + ξ_k(s)
```

where X̃_k(s) is a standardised terrain-distance covariate and ξ_k(s) is a zero-mean Matern GRF with smoothness ν = 1.

**True parameters:**

| Parameter    | Water | Cave | Grass |
|--------------|-------|------|-------|
| β1           | −1.5  | −2.0 | −1.2  |
| ρ (range)    |  2.0  |  1.5 |  2.5  |
| σ            |  0.8  |  1.0 |  0.7  |

**Inference:** INLA-SPDE with Penalised Complexity priors, empirical Bayes integration strategy (`int.strategy = "eb"`), dual-mesh Poisson GLMM reformulation (Simpson et al., 2016).

**Reference coordinates (F_k):** The full coordinate sets used to compute the terrain-distance covariates are tabulated in `outputs/tables/13_appendix_A_reference_coordinates.csv` and in Appendix A of the manuscript.

---

## Key Results (500 replicates)

| Parameter | Guild | 95% CI Coverage | MC SE  |
|-----------|-------|-----------------|--------|
| β1        | Water | 93.0%           | 1.14pp |
| β1        | Cave  | 90.0%           | 1.34pp |
| β1        | Grass | 91.2%           | 1.27pp |
| ρ         | Water | 88.2%           | 1.44pp |
| ρ         | Cave  | 55.8%           | 2.22pp |
| ρ         | Grass | 99.0%           | 0.44pp |
| σ         | Water | 90.2%           | 1.33pp |
| σ         | Cave  | 90.6%           | 1.31pp |
| σ         | Grass | 89.4%           | 1.38pp |

Full results with bias, RMSE, and WAIC discrimination rates are in `outputs/tables/08_simulation_summary_with_mcse.csv`.

---

## Citation

If you use this pipeline, please cite the associated manuscript:

> Irawan, R. (2026). Computational Benchmarking of the INLA-SPDE Approach for Spatial Point Processes: Parameter Recovery and Spatial Confounding. *Journal of Statistical Computation and Simulation* (under revision).

---

## License

This code is made available for academic reproducibility purposes. Please contact the author for any other intended use.
