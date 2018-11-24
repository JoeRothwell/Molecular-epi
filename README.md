# Data analysis for molecular epidemiology projects

The functions of the various scripts stored in this repository are as follows:

`Biocrates_controls_data_prep.R` and `CRC_data_prep.R`: Prepare the EPIC controls dataset and the CRC case-control datasets for analysis, by subsetting the appropriate observations and joining the necessary variables from other datasets.

`CRC_models.R`: Runs conditional logistic models for both score and signature and outputs forest plots. Also meta-analyses data.

`Metabolite_volcano_WCRF.R`: Calculates associations between WCRF scores (high and low categories) and Biocrates metabolites or blood fatty acids by logistic regression.

`Metabolite_signature_WCRF.R`: Finds metabolic signatures of WCRF scores by PLS for Biocrates metabolites or fatty acids. Common metabolites between datasets are first determined and this overlap only used for PLS.

`Plot_signatures.R`: Plots labelled scatter graphs and barcharts to visualise metabolic signatures.

`ukb_analysis.R`: Performs Cox proportional hazards survival models on UK Biobank data to test for associations between diabetes prevalence and CRC. Outputs tidy tables of coefficients and forest plots.

