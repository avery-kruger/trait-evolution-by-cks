library(taxize)
library(ape)

tip.phyla <- data.frame(species=zannetree.trim$tip.label, phylum = NA)
tip.phyla[6814,2] <- "Magnoliopsida" #DB lookup issue for this species
for(i in 6815:nrow(tip.phyla)){#(1:nrow(tip.phyla))[-6814]){
  tip.phyla$phylum[i] <- classification(tip.phyla$species[i], db="gbif", rows = 1)[[1]][3,1]
}

angiosperm.tree <- drop.tip(zannetree.trim, tip.phyla[tip.phyla$phylum != "Magnoliopsida",1])

write.tree(angiosperm.tree, "Zanne.angiosperm.tre")
