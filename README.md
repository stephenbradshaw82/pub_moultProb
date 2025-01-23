# Estimating moult timing in lobsters: A Bayesian Approach

Primary author: Stephen Bradshaw


Contact: stephen.bradshaw@utas.edu.au


Date: 2025-01-23


This respository supports a publication (titled above) to estimate crustacean moult timing using a Hidden Markov model written within a Bayesian framework using STAN (executed using the statistical programming language R).


## Introduction
This repo provides a series of scripts:

(i) 100genDataModel.R 
  - Create simulated data and for varying case studies to test model sensitivity
  - Execute STAN code to estimate moult timing

(ii) 200visualisations.R
  - Sample visualisation code


## Files
```plaintext
Parent Folder
├── 00_src
│   └── functions.R
├── 100genDataModel.R
├── 200visualisations.R
├── models
│   ├── rSTAN_moultProb_np.stan
│   └── rSTAN_moultProb_priors.stan
└── outputs
    ├── ...
    └── plots
```
