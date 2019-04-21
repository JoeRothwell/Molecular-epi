# Data analysis for molecular epidemiology projects

This repository is the codebase for the following projects:

---

#### Metabolic signatures of healthy lifestyle and association with colorectal cancer in EPIC

`Biocrates_controls_data_prep.R`: Plots PCAs of EPIC controls data with adjustment for confounders using the residuals method. 

`CRC_data_prep.R`: Prepares three datasets, the EPIC controls and the two CRC case-control datasets for analysis, by subsetting the appropriate observations and joining the necessary variables from other datasets.

`Biocrates_cmpd_metadata.R`: Extracts coded compound names and full compound names from a Stata .dta file and places them in a data frame.

`Metabolic_signatures.R`: Finds endogenous (Biocrates) and fatty acids signatures of WCRF score using a PLS model and stores these as objects to subsequently use for association with CRC risk. Common metabolites between datasets are first determined and this overlap only used for PLS.

`CRC_models_calc_score.R` and `CRC_models_signature.R`: Run Cox or conditional logistic models for score and signature respectively. Create forest plots and meta-analyse data using `metafor`.

`Metabolite_volcano_WCRF.R`: Calculates associations between WCRF scores (high and low categories) and Biocrates metabolites or blood fatty acids by logistic regression.

`Plot_signatures.R`: Plots labelled scatter graphs and barcharts to visualise metabolic signatures.

`baseline_characteristics.R`: Generates Table 1 automatically from raw data using the `qwraps2` package.

`GRS_biocrates.R` and `GRS_fatty.acids.R`: Calculate associations between polygenic risk scores of CRC and biocrates metabolites or fatty acids.

---

#### Association of type 2 diabetes and colorectal cancer in the UK Biobank cohort

`ukb_analysis.R`: Performs Cox proportional hazards survival models on UK Biobank data to test for associations between diabetes prevalence and CRC. Outputs tidy tables of coefficients and forest plots.

---

#### A metabolomics study on coffee consumption and liver cancer in the ATBC study

`ATBC_case_control.R` and `ATBC_metabolomics.R`: subsetting and feature filtering of raw metabolomics data for the ATBC coffee study.

---

#### Association of serum metabolites and dietary and lifestyle factors in the Lifepath study

`Lifepath_exploratory.R`: Exploratory analysis of Lifepath NMR metabolomics data.

`baseline_Lifepath.Rmd`: Generates a baseline characteristics file for the study from subject metadata.

