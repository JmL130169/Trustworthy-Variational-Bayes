# Bayesian Mixture Linear Regression (BMLR)

This folder contains the implementation of Trustworthy Variational Bayes (TVB) for Bayesian Mixture Linear Regression models as described in "Bend to Mend: Toward Trustworthy Variational Bayes with Valid Uncertainty Quantification" by Jiaming Liu and Meng Li.

## Overview

Bayesian Mixture Linear Regression models are used when the data is believed to come from multiple subpopulations, each following its own linear regression model. This implementation provides:

1. Standard Variational Bayes (VB) for BMLR
2. Fractional VB for uncertainty quantification calibration
3. Two methods for TVB calibration: sequential update and grid-search approach

## Key Functions

- `vb_bmlr`: Main variational Bayes function for BMLR models
- `update_omega_beta`: Sequential TVB update function
- `boot_tvb`: Grid search dictionary building function

## Simulation Parameters
The simulation studies in the paper used the following configuration:
```bash
--name_batch sim601 --seed 1-500 --cluster 2 --dataset 500 --dim 2 --SNR 0.1 --mixing 0.65
