# Data analysis for molecular epidemiology projects

The functions of the various scripts stored in this repository are as follows:

`Biocrates_controls_data_prep.R` and `CRC_data_prep.R`: Prepare the EPIC controls dataset and the CRC case-control datasets for analysis, by subsetting the appropriate observations and joining the necessary variables from other datasets.

`Metabolite_volcano_WCRF.R`: Calculates associations between WCRF scores (high and low categories) and Biocrates metabolites or blood fatty acids by logistic regression.

`Metabolite_signature_WCRF.R`: Finds metabolic signatures of WCRF scores by PLS for Biocrates metabolites or fatty acids. Models CRC case-control status by metabolic signature or calculate score (conditional logistic regression) and plots the output as a forest plot. Also meta-analyses large and small CRC case-controls.

`Plot_signatures`: Plots labelled scatter graphs and barcharts to visualise metabolic signatures.

`ukb_analysis.R`: Performs Cox proportional hazards survival models on UK Biobank data to test for associations between diabetes prevalence and CRC.

