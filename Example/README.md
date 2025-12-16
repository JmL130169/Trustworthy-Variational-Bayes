# Example: Old Faithful Dataset Analysis

This folder contains example applications of Trustworthy Variational Bayes (TVB) to the Old Faithful geyser dataset using Gaussian Mixture Models (GMM).

## Overview

The Old Faithful dataset, included in base R, provides observations of eruption durations and waiting times for the Old Faithful geyser in Yellowstone National Park. These eruptions naturally form two clusters, making it an excellent candidate for applying Gaussian Mixture Models. This example demonstrates:

1. Application of standard VB to a real dataset
2. Implementation of both sequential and grid-search TVB approaches

## Files

- `TVB_support_gmm.R`: Support functions for fractional VB estimation and TVB calibration
- `realdata_GMM.R`: Main script for analyzing the Old Faithful dataset
- `gmm_realdata.RData`: Estimation results
- `realdata_plot.R`: Uses "gmm_realdata.RData" as dictionary to reproduce the plot in Section 7 of the paper.

## Key Functions

- `boot_gmm`: Function to estimate all Bootstrap datasets for a specific $\omega$
- `omega_bootstrap`: Function for building dictionary over all $\omega$'s, used with `boot_gmm`
- `update_omega_beta`: Sequential update function for the sum of mean vectors
- `update_omega_pi`: Sequential update function for mixing probabilities

## Main Analysis Script: `realdata_GMM.R`

This script contains the complete analysis workflow for applying TVB to the Old Faithful dataset using Gaussian Mixture Models.

### Workflow:

1. **Data Loading and Preparation**:
  - Loads the built-in `faithful` dataset in R
  
2. **Dictionary Building for Grid-Search TVB**:
  - Creates a fine grid of ω values from 0.001 to 1
  - Calls `omega_bootstrap` function to build the calibration dictionary
  - For each $\omega$ value, computes VB estimates and bootstrap samples
  - Stores all results in a comprehensive dictionary structure

3. **Grid-Search TVB Calibration**:
  - Uses the dictionary to identify optimal $\omega$ values for:
    - Mixing probability ($\pi$)
    - Sum of component means in the first cluster ($\beta_{1}$)
  - Applies the optimal $\omega$ values to construct calibrated credible intervals

4. **Sequential TVB Implementation**:
  - Runs `update_omega_pi` for mixing probability calibration
  - Runs `update_omega_beta` for sum of means calibration

5. **Standard VB Baseline**:
  - Computes standard VB estimates and credible intervals
  - Uses these as a baseline for comparison

6. **Result Storage**:
  - Saves all analysis results to "gmm_realdata.RData"

This comprehensive script serves as a practical example of how to apply TVB methods to real-world data, demonstrating both the sequential and grid-search approaches to uncertainty quantification calibration.

## References

For more details on the methodology, please refer to the main paper "Bend to Mend: Toward Trustworthy Variational Bayes with Valid Uncertainty Quantification" by Jiaming Liu and Meng Li. The Old Faithful dataset is described in:
- Azzalini, A., & Bowman, A. W. (1990). A look at some data on the Old Faithful geyser. Journal of the Royal Statistical Society. Series C (Applied Statistics), 39(3), 357-365.
