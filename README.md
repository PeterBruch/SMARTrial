
# SMARTrial

This repository contains data and executable scripts to reproduce the figures and analyses presented in the paper: 

_**Ex vivo drug response profiling for response and outcome prediction in hematologic malignancies: The prospective non-interventional SMARTrial**_

_Nora Liebers\*, Peter-Martin Bruch\*, Tobias Terzer, Miguel Hernandez, Carolin Kolb, Mareike Knoll, Angela Lenze, Donnacha Fitzgerald, Carsten Mueller-Tidow, Wolfgang Huber, Junyan Lu, Axel Benner, Thorsten Zenz, Sascha Dietrich_

\*These authors contributed equally to this work. 

All relevant data is stored in `/data/SMARTrial_data.RData` and is automatically loaded in the analysis script. The analysis script can be found under `vignettes/src/SMARTrial_Analysis.Rmd` and the associated .html file (showing the analysis in full) in `inst/doc/SMARTrial_Analysis.html`. To run the entire analysis, clone the repository and run the script `SMARTrial_Analysis.Rmd`. Make sure you have updated R to at least version 4.1.0 and have installed all libraries in `SMARTrial_Analysis.Rmd` beforehand. 

If you use this analysis in published research, please cite the paper. Please refer to the manuscript for more details on experimental methods and analysis. 

_***Installation guide***_
To run the entire analysis, clone the repository and run the script `SMARTrial_Analysis.Rmd`. 

_***System requirements***_
To run the analysis, a R (at least version 4.1.0) and all dependency libraries specified under `inst/doc/SMARTrial_Analysis.html#16_16_Session_info` are required. 

_***Output***_
The analysis script takes roughly 30 minutes to run, mostly due to the elastic net analysis (Section 8) and Drug~Drug correlations (Section 12). The expected output is shown in `inst/doc/SMARTrial_Analysis.html`.  