# Trustworthy Variational Bayes (TVB)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the implementation of "Bend to Mend: Toward Trustworthy Variational Bayes with Valid Uncertainty Quantification" by Jiaming Liu and Meng Li (Rice University) https://arxiv.org/abs/2512.22655.

## Overview

Variational Bayes (VB) is a computationally efficient method for approximating posterior distributions in Bayesian inference. However, VB often provides unreliable uncertainty quantification, particularly in the form of credible intervals that fail to achieve nominal coverage. This repository implements Trustworthy Variational Bayes (TVB), a novel approach to recalibrate VB procedures and provide valid uncertainty quantification.

### Key Features

- **Bend-to-mend strategy**: Intentionally misspecify the likelihood (bend) to correct VB's flawed uncertainty quantification (mend)
- **Two calibration methods**: 
  - Sequential update TVB: Iteratively adjusts the calibration parameter
  - Grid-search TVB: More efficient approach utilizing a pre-computed dictionary of calibration parameters
- **Theoretical guarantees**: Calibrated credible intervals achieve correct frequentist coverage
- **Applications**: Implemented for Gaussian mixture models (GMM) and Bayesian mixture linear regression (BMLR)

## Repository Structure

| File/Directory | Description |
|----------------|-------------|
| **GMM/** | Gaussian Mixture Model |
| &nbsp;&nbsp;&nbsp;&nbsp;tvb_gmm.R | Functions and simulations of GMM |
| &nbsp;&nbsp;&nbsp;&nbsp;README.md | Detail about GMM implementation |
| **BMLR/** | Bayesian Mixture Linear Regression |
| &nbsp;&nbsp;&nbsp;&nbsp;tvb_bmlr.R | Functions and simulations of BMLR |
| &nbsp;&nbsp;&nbsp;&nbsp;README.md | Detail about BMLR implementation |
| **Example/** | Example applications on real data (GMM) |
| &nbsp;&nbsp;&nbsp;&nbsp;TVB_support_gmm.R | Support functions of using GMM on real data |
| &nbsp;&nbsp;&nbsp;&nbsp;realdata_GMM.R | Implementation of GMM on `old faithful` dataset |
| &nbsp;&nbsp;&nbsp;&nbsp;gmm_realdata.RData | Results of real data application |
| &nbsp;&nbsp;&nbsp;&nbsp;README.md | Detail about real data example implementation |
| README.md | This file |

## Citation
```bibtex
@article{liu2025bend,
  title={Bend to Mend: Toward Trustworthy Variational Bayes with Valid Uncertainty Quantification},
  author={Liu, Jiaming and Li, Meng},
  journal={arXiv preprint arXiv:2512.22655},
  year={2025}
}
```

## License
This project is licensed under the MIT License - see the LICENSE file for details.


