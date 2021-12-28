# GribenskiAlps2022

This repository contains the Matlab code needed to reproduce the results in the following manuscript: 

Gribenski, N., Valla, P.G. Tremblay, M.M., Guralnik, B., Balco, G., and Shuster, D.L. Cosmogenic 3He paleothermometry on post-LGM glacial bedrock within the central European Alps. To be submitted to the journal Geochronology. 

There are two different sets of code associated with this manuscript. The first set of code is located in the folder "MDDDiffusionKinetics". The code in this folder can be used to model noble gas diffusion kinetics from laboratory step heating experiments using a multiple diffusion domain (MDD) model. The user should edit inputs and run the code "SetupMDD.m". An example dataset for implementing the MDD model called "MBTP9-MDD.txt" is provided in folder "Inputs". 

The second set of code is located in the folder "ProductionDiffusionForwardSimulation". The code in this folder can be used to run forward simulations of cosmogenic noble gas production and diffusion for different exposure and thermal history scenarios. The user should edit inputs and run the code "wrap_PD_MDD_LGMScenario_MBTP1_fig3". An example dataset for implementing the forward model called is provided in folder "Data", and the user can reproduce Figure 3 using this example dataset and the example thermal history provided in the main folder, titled "FixedLGMmin15_MBTP1.txt". 

Any questions about these Matlab codes should be directed to Marissa Tremblay, tremblam@purdue.edu. 

[![DOI](https://zenodo.org/badge/442584835.svg)](https://zenodo.org/badge/latestdoi/442584835)
