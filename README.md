# Learning Virulence–Transmission Relationships Using Causal Inference

This repository contains MATLAB code and data supporting the analyses presented in the article:

**"Learning Virulence–Transmission Relationships Using Causal Inference"**

## Overview

The MATLAB scripts provided here facilitate the analysis and visualization of epidemiological data, focusing on the relationship between virulence ($\mu$) and transmission ($\beta$). The repository includes methods for estimating causal relationships, computing correlation dimensions, and the computation and visualization of invariant densities.

## Features

The scripts in this repository enable:

- Processing COVID-19 epidemiological datasets.
- Estimation of virulence ($\mu$) and transmission ($\beta$) dynamics.
- Calculation of correlation dimensions and geometric causal measures.
- Estimation of invariant densities.
- Generation of figures for time series, density distributions, and causal relationships.

### Key Scripts

- **COVIDRunnerV2.m / EuropeCOVIDRunner.m / USACOVIDRunner.m / AsiaCOVIDRunner.m** – COVID-19 data preprocessing, visualization, and computation of geometric causal measures.  
- **estimateInvariantDensity.m** – Computes invariant densities for virulence or transmission variables.  
- **Runner2.m / RunnerCd.m / densityPlots.m** – Synthetic data generation, computation of geometric causal measures, and visualization of invariant densities.

## Usage

1. Load the relevant `.mat` data file corresponding to your region or dataset.  
2. Run the appropriate MATLAB script to compute measures or generate plots.  
3. Outputs include figures and CSV files containing geometric causal measures.

## Citation

If you use this code for research purposes, please cite the original article:

*"Learning Virulence–Transmission Relationships Using Causal Inference"*
