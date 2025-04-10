# Mosquito-suppression-trials
This repository contains code for generating simulated data from cluster randomised control trials (CRCTs) aiming to detect suppression of malaria vector populations. The code implements the methodology described in:

P.A. Hancock, T. J. Hui, P. S. Epopa, A. Milogo, A. R. McKemey, F. A. Yao,  A. Diabaté, A. Burt, "Requirements for designing cluster randomised control trials to detect suppression of malaria vector population densities" 

Individual directories contain code written in R and Stata, and the [associated vignette](https://pahanc.github.io/Mosquito-suppression-trials/) describes how the code is implemented to produce the results presented in the study. 

System requirements: This code is written in R and Stata software, which runs on a wide variety of UNIX platforms, Windows and Mac OS. The R software is quick to install. Stata can be purchased from https://www.stata.com/order/.

Software requirements: R, using the following packages: R-INLA, data.table, Matrix, scales. This code has been tested on R version 4.4.2 with the package R-INLA version 24.12.11, data.table version 1.16.4, Matrix version 1.7-1, scales version 1.3.0. Stata code has been tested on Stata version 18.5/
