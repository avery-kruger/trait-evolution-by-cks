library(taxize)
library(ape)

zannetree <- read.tree('Vascular_Plants_rooted.dated.tre')

zannetree <- read.tree("Zanne.angiosperm.tre")
zannefreeze <- read.csv("MinimumFreezingExposure.csv")
zannetree$tip.label <- trimws(tolower(zannetree$tip.label),"both")
zannefreeze$species <- trimws(gsub(" ","_",tolower(zannefreeze$species)),"both")
nofreezenames <- which(is.na(match(zannetree$tip.label,zannefreeze$species)))
myfreeze <- zannefreeze[which(!is.na(match(zannefreeze$species,zannetree$tip.label))),]
zannetree.trim <- drop.tip(zannetree,nofreezenames)

tip.phyla <- data.frame(species=zannetree.trim$tip.label, phylum = NA)
#Uncomment to find phyla for each species. Has some database lookup issues.
# for(i in 1:nrow(tip.phyla)){
#   myname <- tip.phyla$species[i]
#   titlename <- paste0(toupper(substr(myname,1,1)), substr(myname,2, nchar(myname)))
#   titlename <- gsub("_"," ", titlename)
#   tip.phyla$phylum[i] <- classification(titlename, db="gbif", rows = 1)[[1]][3,1]
# }
tip.phyla$phylum[3579:13381] <- "Magnoliopsida"

angiosperm.tree <- drop.tip(zannetree.trim, tip.phyla[tip.phyla$phylum != "Magnoliopsida",1])

write.tree(angiosperm.tree, "Zanne.angiosperm.tre")
