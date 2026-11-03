# LETR code

MATLAB R2019b or later (`clim` in the figures needs R2022a or later) with the Statistics and Machine Learning Toolbox. Run the scripts from this folder. The SARS-CoV-2 data are in `COVID/`. Each runner writes its outputs to its own `Results...` folder.

## Runners

| File | Purpose |
|---|---|
| `runCase1a.m` | Runs Case 1(a), direct coupling, and writes the GeoC table, figures and the alternative-methods table. |
| `runCase1b.m` | Runs Case 1(b), variants (i) to (iii) with a shared latent driver, and writes the GeoC tables and figures. |
| `runCase1c.m` | Runs Case 1(c), the myxomatosis chain, and writes the GeoC table, invariant density and overview figure. |
| `runCOVID.m` | Runs the SARS-CoV-2 analysis for the global, Asia, Europe and USA records over the three periods. |
| `runSensitivity.m` | Runs the one-at-a-time sensitivity of GeoC to Cmax and to the local-slope half-window. |
| `runSampleSize.m` | Runs the sensitivity of GeoC to the number of transitions N and draws the Case 1(a) local-slope plots for several N. |

## Core functions

| File | Purpose |
|---|---|
| `estimateD2.m` | Estimates the correlation dimension from a plateau of the local-slope profile. |
| `conditionalD2.m` | Computes the conditional dimension D2(X given Y) as D2([Y,X]) minus D2(Y). |
| `estimateGeoC.m` | Computes both directional GeoC values and returns a table with the details for plotting. |
| `estimateGeoCTI.m` | Computes a compact GeoC that returns numbers only and is used by the sensitivity runner. |
| `estimateInvariantDensity2D.m` | Builds the boundary-reflected transfer operator and returns its invariant density. |

## Data generation and preparation

| File | Purpose |
|---|---|
| `generateCase1Series.m` | Generates the deterministic benchmark records of Case 1. |
| `dwyerBeta.m` | Computes the transmission of a strain from the within-host model of Dwyer et al. |
| `defineCOVIDRegionConfigs.m` | Lists the regions, data files and period ranges of the SARS-CoV-2 analysis. |
| `preprocessCOVIDPeriod.m` | Selects a period, removes non-finite pairs and min-max scales Beta within the period. |
| `buildSensitivitySystems.m` | Collects the benchmark and SARS-CoV-2 records analyzed by the sensitivity runner. |

## SARS-CoV-2 workflow

| File | Purpose |
|---|---|
| `analyzeCOVIDRegionLETR.m` | Runs GeoC and the invariant density for one region over the three periods. |
| `saveCOVIDRegionResults.m` | Saves the result of one region as a MAT file and one CSV table per period. |
| `makeCOVIDRegionFigures.m` | Saves the local-slope figure and the invariant-density figure of one region. |
| `makeCOVIDRegionalMainFigure.m` | Draws the regional synthesis figure for Asia, Europe and USA. |

## Figures

| File | Purpose |
|---|---|
| `drawLocalSlopePanel.m` | Draws one local-slope panel with two D2 profiles, scaling windows and a GeoC box. |
| `plotLocalSlopeDiagnostics2x2.m` | Draws the 2-by-2 local-slope figure for the four GeoC terms of one record. |
| `plotLocalSlopeDiagnostics3Cases.m` | Draws the 3-by-4 local-slope figure with one row per case. |
| `plotLocalSlopeAcrossN.m` | Draws the local-slope comparison of D2(beta_n) and D2([beta_n,beta_{n+1}]) for several values of N. |
| `plotInvariantDensity3Periods.m` | Draws the 3-by-3 invariant-density figure of one region. |
| `addPanelLabel.m` | Adds a bold panel label such as (a) to an axes. |
| `batlow.m` | Returns the perceptually uniform batlow colormap. |

## Alternative methods

| File | Purpose |
|---|---|
| `transferEntropyKSG.m` | Estimates one-step transfer entropy with a k-nearest-neighbor estimator and a permutation test. |
| `ccmSkill.m` | Computes convergent cross mapping skill by simplex projection. |
| `grangerLinearF.m` | Runs the classical linear Granger F test. |
| `makeAltMethodsLatexTable.m` | Writes a LaTeX table comparing GeoC with transfer entropy, CCM and Granger. |

The header of every file lists its calls, callers, inputs and outputs.
