# glioma-Tcell-MDSC

<a href="https://github.com/stepien-lab/glioma-Tcell-MDSC/"><img src="https://img.shields.io/badge/GitHub-stepien--lab%2Fglioma--Tcell--MDSC-blue" /></a> <a href="https://doi.org/10.1101/2023.05.15.540846"><img src="https://img.shields.io/badge/bioRxiv-2023.05.15.540846-orange" /></a> <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" /></a>

The code contained in the glioma-Tcell-MDSC project was developed to numerically simulate the interactions of glioma cells, T cells, and myeloid-derived suppressor cells (MDSCs) in a model of brain cancer. It is described in:
>[Hannah G. Anderson](https://github.com/HannahGrace314), [Gregory P. Takacs](https://pharmacology.med.ufl.edu/profile/takacs-gregory/), [Duane C. Harris](https://search.asu.edu/profile/2524814), [Yang Kuang](https://math.la.asu.edu/~kuang/), [Jeffrey K. Harrison](https://pharmacology.med.ufl.edu/profile/harrison-jeffrey/), and [Tracy L. Stepien](https://github.com/tstepien/), Global stability and parameter analysis reinforce therapeutic targets of PD-L1-PD-1 and MDSCs for glioblastoma, Submitted to _Journal of Mathematical Biology_, bioRxiv: [2023.05.15.540846](https://doi.org/10.1101/2023.05.15.540846).

## Programs
+ [GBMFunc.m](GBMFunc.m): model equations
+ [ABCrejection.m](ABCrejection.m): run this program to obtain parameter samples and error values according to data
+ [ABCfigures.m](ABCfigures.m): accepts parameter values based on specified error threshold to finish ABC rejection and obtain scatter/contour plots
+ [GBMnumsim.m](GBMnumsim.m): plots numerical simulations of model
+ [ABCdistributions.m](ABCdistributions.m): determines distributions of ABC histograms
+ [GBMbifurcationanalysis.m](GBMbifurcationanalysis.m): produces bifurcation analysis on parameters

![image](https://user-images.githubusercontent.com/89090482/209017684-ac768527-f079-4604-a4ca-e719dde711b5.png)


Code for the extended Fourier Analysis Sensitivity Test (eFAST) method was developed by the Krishner lab at Univeristy of Michigan and can be obtained at:

Kirschner D (2008) Uncertainty and sensitivity functions and implementation.
http://malthus.micro.med.umich.edu/lab/usadata/


## Lead Developer
The lead developer of this code is [Hannah G. Anderson](https://github.com/HannahGrace314).

## Licensing
Copyright 2022-2023 [Tracy Stepien's Research Lab Group](https://github.com/stepien-lab/). This is free software made available under the MIT License. For details see the [LICENSE](LICENSE) file.
