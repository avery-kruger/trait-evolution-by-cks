library(taxize)
library(ape)

zannetree <- read.tree("Zanne.angiosperm.tre")
zannefreeze <- read.csv("MinimumFreezingExposure.csv")
zannetree$tip.label <- trimws(tolower(zannetree$tip.label),"both")
zannefreeze$species <- trimws(gsub(" ","_",tolower(zannefreeze$species)),"both")
nofreezenames <- which(is.na(match(zannetree$tip.label,zannefreeze$species)))
myfreeze <- zannefreeze[which(!is.na(match(zannefreeze$species,zannetree$tip.label))),]
zannetree.trim <- drop.tip(zannetree,nofreezenames)

tip.phyla <- data.frame(species=zannetree.trim$tip.label, phylum = NA)
tip.phyla[6814,2] <- "Magnoliopsida" #DB lookup issue for this species
for(i in 6815:nrow(tip.phyla)){
  tip.phyla$phylum[i] <- classification(tip.phyla$species[i], db="gbif", rows = 1)[[1]][3,1]
}

angiosperm.tree <- drop.tip(zannetree.trim, tip.phyla[tip.phyla$phylum != "Magnoliopsida",1])

write.tree(angiosperm.tree, "Zanne.angiosperm.tre")
