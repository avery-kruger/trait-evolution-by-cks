README for Inferring the Evolutionary Model of Community-Structuring Traits with Convolutional Kitchen Sinks

Install kitchen R package:
library(devtools)
install_github("avery-kruger/kitchen")

Simulated data are located in data.zip

Simulated Data
To reproduce analysis of convolutional kitchen sinks (CKS) performance for predicting traits governing community assembly by filtering:
(Optional: Run LimSim_simulations.R to simulate data._
Run LimSim_cks.R to reproduce analysis. 

To reproduce analysis of convolutional kitchen sinks (CKS) performance for predicting traits governing community assembly by limiting similarity:
(Optional: Run LimSim_simulations.R to simulate data.)
Run LimSim_cks.R to reproduce analysis.

To reproduce empirical analysis:
Download the Zanne et al. dataset from https://doi.org/10.5061/dryad.63q27. 
Move the files MinimumFreezingExposure.csv (in the climate folder) and Vascular_Plants_rooted.dated.tre (in the Phylogenetic Resources folder) to your working directory.
Run Zanne_tree.R to trim the Zanne tree to only contain angiosperms represented in MinimumFreezingExposure.csv.(Optional: Run Zanne_simulations.R to simulate communities based on the Zanne phylogeny.)
Run Zanne_cks.R to conduct CKS and MPD analysis on the simulated and empirical freeze-tolerant communities.

