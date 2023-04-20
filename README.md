
# SMARTrial

This repository contains data and executable scripts to reproduce the figures and analyses presented in the paper: 

_**Ex vivo drug response profiling for response and outcome prediction in hematologic malignancies: The prospective non-interventional SMARTrial**_

_Nora Liebers\*, Peter-Martin Bruch\*, Tobias Terzer, Miguel Hernandez-Hernandez, Nagarajan Paramasivam, Donnacha Fitzgerald, Heidi Altmann, Tobias Roider, Carolin Kolb, Mareike Knoll, Angela Lenze, Uwe Platzbecker, Christoph Röllig, Claudia Baldus, Hubert Serve, Martin Bornhäuser, Daniel Hübschmann, Carsten Müller-Tidow, Friedrich Stölzel, Wolfgang Huber, Axel Benner, Thorsten Zenz, Junyan Lu, Sascha Dietrich_

\*These authors contributed equally to this work. 

All data of the initial SMARTrial cohort and the AML validation cohort is stored in `/data/SMARTrial_data.RData` and `/data/Validation_data.RData`, respectively. The analysis scripts can be found under `vignettes/src/SMARTrial_Analysis.Rmd` and `vignettes/src/Validation_Analysis.Rmd` along with the associated .html files (showing the analysis in full). 

If you use this analysis in published research, please cite the paper. Please refer to the manuscript for more details on experimental methods and analysis. 

_***Installation guide***_
To run the entire analysis, clone the repository and run the script `SMARTrial_Analysis.Rmd` and `Validation_Analysis.Rmd`. 

_***System requirements***_
To run the analysis, R (at least version 4.1.0) and all dependency libraries are required. 

_***Output***_
The SMARTrial analysis script takes roughly 30 minutes to run, mostly due to the elastic net analysis (Section 8) and Drug~Drug correlations (Section 12). The expected output is shown in `vignettes/src/SMARTrial_Analysis.html` and `vignettes/src/Validation_Analysis.html`.  