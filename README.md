# Long-term impact of pneumococcal conjugate vaccines on invasive pneumococcal disease incidence among all ages from national, active, laboratory-based surveillance, South Africa, 2005-2019
This repository houses data used for the von Gottberg et al. 2024 manuscript.
Created by: Jackie Kleynhans, National Institute for Communicable Diseases.
We evaluated the ongoing effects (direct and indirect) of prevention of pneumococcal conjugate vaccine (PCV) serotypes, replacement disease, and effects on antimicrobial resistance, over 15 years of sustained surveillance and before the COVID-19 pandemic. In South Africa, PCV7 was introduced in 2009, and PCV13 in 2011, in a two plus one schedule. We conducted national, active, laboratory-based surveillance for invasive pneumococcal disease (IPD) among all ages, including isolate serotyping and susceptibility testing through GERMS-SA. We fitted individual linear regression models to imputed IPD case counts from 2005 through 2019 by serotype and age to estimate and compare expected and actual IPD cases in 2019.

The analytical framework was adapted from existing code for the evaluation of PCV impact: Weinberger DM. No Title. isppd_workshop2022. GitHub. Published 2022. https://github.com/DanWeinberger/isppd_workshop2022/tree/main.

All analyses were performed in R studio version 23.12.0.

## Included in this repository:
### Datasets: 
Annual and quarterly imputed IPD counts for 2005-2019 (Annual imputed IPD counts 2024-02-20.csv and Quarterly imputed IPD counts 2024-02-20.csv)
### Scripts: 
Analysis: PCV Impact ITS 2005-2019 2024-02-20.R
Function for ITS: functions its 2024-02-20.R (sourced from analysis script)
Plotting overall incidences: PCV Impact incidence plots 2005-2019 2024-02-20.R
### Results:
Results are ordered in the ITS results and plots folder, based on aggregation level (annual/quarterly), distribution (negative binomial and Poisson) and if a cosine term was used for seasonality. The main analysis presented in the manuscript is in the Annual_Negative_Binomial folder. All other model results are provided as supplementary information. 
