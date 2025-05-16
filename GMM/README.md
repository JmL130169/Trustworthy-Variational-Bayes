# Gaussian Mixture Model (GMM)

This folder contains the implementation of Trustworthy Variational Bayes (TVB) for Gaussian Mixture Models as described in "Bend to Mend: Toward Trustworthy Variational Bayes with Valid Uncertainty Quantification" by Jiaming Liu and Meng Li.

## Overview

Gaussian Mixture Models are widely used for clustering and density estimation when data is believed to come from multiple subpopulations, each following a Gaussian distribution. This implementation provides:

1. Standard Variational Bayes (VB) for GMM
2. Fractional VB for uncertainty quantification calibration
3. Two methods for TVB calibration: sequential update and grid-search approach

## Key Functions

- `vb_gmm`: Main variational Bayes function for GMM models
- `update_omega_pi`: Sequential TVB update function for mixing probability
- `boot_tvb`: Grid search dictionary building function

## Simulation Parameters
The simulation studies in the paper used the following configurations:
For `n = 1000`
```bash
--name_batch sim306 --seed $SLURM_ARRAY_TASK_ID --cluster 2 --dim 2 --center 2 --sample 1000 --prob 0.65
```
where `$SLURM_ARRAY_TASK_ID` takes value `1 - 500`. For `n = 1500, 2000, 2500, 3000`, change `--name_batch` to `sim308, sim310, sim312, sim314`, and `--sample`, respectively.
