library(geiger)
library(picante)
library(phytools)
library(PhyloMeasures)
#library(ecolottery)
#library(vegan)
library(reshape2)
library(parallel)
source("sim_MPD_functions.R")

zannetree <- read.tree("Zanne.angiosperm.tre")
zannefreeze <- read.csv("MinimumFreezingExposure.csv")
zannetree$tip.label <- trimws(tolower(zannetree$tip.label),"both")
zannefreeze$species <- trimws(gsub(" ","_",tolower(zannefreeze$species)),"both")
nofreezenames <- which(is.na(match(zannetree$tip.label,zannefreeze$species)))
myfreeze <- zannefreeze[which(!is.na(match(zannefreeze$species,zannetree$tip.label))),]
zannetree.trim <- drop.tip(zannetree,nofreezenames)

commsize <- sum(myfreeze$Freeze.tmin.lo == "FreezingExposed")

#Broken into 200 25-sim chunks for computational reasons
for(i in 1:200){
  newsim <- community_sim(zannetree.trim,
                  reps=25,
                  comm_size = commsize,
                  param_list = list(n=1,mean=0,sd=0.08),
                  param_bounds = c(-.2,.2),
                  optimum_species = which(myfreeze$Freeze.tmin.lo == "FreezingExposed")
                  )
  saveRDS(newsim, file = paste0("angiocomm",i,".rds"))
  print(paste0(i,"/",200, " Community Simulated"))
  mpd <- MPDapply(newsim$comm, tree=zannetree.trim)
  mpd$param <- newsim$param
  saveRDS(mpd, file = paste0("angiosim",i,".rds"))
  print(paste0(i,"/",200, " MPD Calculated"))
}


