library(geiger)
library(picante)
library(phytools)
library(PhyloMeasures)
library(ecolottery)
#library(vegan)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(visreg)
library(ggtree)
library(kitchen)
#load_all("/Users/macuser/kitchen") #replace with regular package load when possible
source("sim_MPD_functions.R")

alldata <- data.frame(true = NA, predicted = NA)
allR2 <- data.frame(ntaxa = c(48, 2^(6:10)), predicted = NA)
for(ntaxa in c(48, 2^(6:10))){
  filename <- paste0(ntaxa,"treesim")
  allsims <- list.files(path="data",pattern=paste0(filename,"_mpd"), full.names = T) %>%
    map_dfr(readRDS)

  allcomm <- list.files(path="data",pattern=paste0(filename,"_comm"), full.names = T) %>%
    map_dfr(readRDS)

  saveRDS(allcomm, file = paste0("allcomm_", ntaxa))
  saveRDS(allsims, file = paste0("allmpd_", ntaxa))


  trainingdata <- allcomm[1:1200,]
  testdata <- allcomm[1201:2000,]

  Y <- trainingdata$param
  X <- trainingdata$comm #[,-which(names(trainingdata) == "param")]
  valY <- testdata$param
  valX <- testdata$comm #[,-which(names(trainingdata) == "param")]

  #Choose superparameters to sweep across
  featuresweep <- 2^(4:9)
  windowsweep <- c(5, 10, 20, 40, 80, 128)
  windowsweep <- windowsweep[windowsweep <= ntaxa]
  if(max(windowsweep) < ntaxa){windowsweep <- c(windowsweep, ntaxa)}

#Conduct preliminary sweep
mysweep <- kitchen_sweep(X,Y,
                         valX,valY,
                         featuresweep,windowsweep,
                         verbose = T, show.plot = T)

#bestR2
best.index <- which(mysweep == max(mysweep), arr.ind = T)
featuresweep[best.index[1]]
windowsweep[best.index[2]]

###Use the best performing model to predict values for the validation data
bestmodelpredictions <- kitchen_prediction(X,
                                           Y,
                                           valX,
                                           featuresweep[best.index[1]],
                                           windowsweep[best.index[2]])

mydata <- data.frame(true = valY, predicted = bestmodelpredictions[[1]][[1]])

alldata <- rbind(alldata, mydata)

###Examine how predictions line up with true data.
ggplot(mydata,
       aes(y=true,x=predicted)) +
  theme_bw() +
  geom_point(aes(alpha = 0.2)) +
  #geom_errorbar(aes(xmin=CIlow95, xmax=CIhigh95))+
  geom_abline(slope=1, intercept=0, color = "blue") +
  xlab("Predicted Parameter") +
  ylab("True Parameter") +
  theme(legend.position = "none")
ggsave(paste0("RandTree_",ntaxa,"SpTree.png"))

allR2[allR2$ntaxa == ntaxa,2] <- summary(lm(mydata$true ~ mydata$predicted))$adj.r.squared
}


#Now compare with MPD
alldata.mpd <- data.frame(true = NA, predicted = NA)
allR2.mpd <- data.frame(ntaxa = c(48, 2^(6:10)), predicted = NA)

for(ntaxa in c(48, 2^(6:10))){
  filename <- paste0(ntaxa,"treesim")
  allmpd <- list.files(path="data",pattern=paste0(filename,"_mpd"), full.names = T) %>%
    map_dfr(readRDS)

  trainingdata.mpd <- allmpd[1:1200,]
  testdata.mpd <- allmpd[1201:2000,]

  Y <- trainingdata.mpd$param
  X <- as.data.frame(trainingdata.mpd$mpd)
  valY <- testdata.mpd$param
  valX <- as.data.frame(testdata.mpd$mpd)

  mymodel <- lm(Y ~ ., data = X)
  mypredictions <- predict(mymodel, as.data.frame(valX))

  mydata <- data.frame(true = valY, predicted = mypredictions)

  alldata.mpd <- rbind(alldata.mpd, mydata)

  ###Examine how predictions line up with true data.
  ggplot(mydata,
         aes(y=true,x=predicted)) +
    theme_bw() +
    geom_point(aes(alpha = 0.2)) +
    #geom_errorbar(aes(xmin=CIlow95, xmax=CIhigh95))+
    geom_abline(slope=1, intercept=0, color = "blue") +
    xlab("Predicted Parameter") +
    ylab("True Parameter") +
    theme(legend.position = "none")
  ggsave(paste0("RandTreeMPD_",ntaxa,"SpTree.png"))

  allR2.mpd[allR2.mpd$ntaxa == ntaxa,2] <- summary(lm(mydata$true ~ mydata$predicted))$adj.r.squared
}

alldata <- alldata[2:4801,]
alldata.mpd <- alldata.mpd[2:4801,]
alldata$ntaxa <- sort(rep(c(48, 2^(6:10)), 800))
alldata.mpd$ntaxa <- sort(rep(c(48, 2^(6:10)), 800))
alldata$type <- "CKS"
alldata.mpd$type <- "MPD"

alldata.both <- rbind(alldata, alldata.mpd)

myR2 <- data.frame(ntaxa = rep(c(48,2^(6:10)),2),
                   type = c(rep("CKS", 6),rep("MPD", 6)),
                   Rsq = NA)
myR2$Rsq[1:6] <- c("0.64","0.66","0.79","0.74","0.76","0.81")
myR2$Rsq[7:12] <- round(allR2.mpd$predicted, 2)
myR2$Rsqlab <- sapply(myR2$Rsq, function(x){paste("bold(","{R^2}==",x,")")})

ggplot(alldata.both,
       aes(y=true,x=predicted, group=ntaxa)) +
  theme_bw() +
  geom_point(aes(alpha = 0.1)) +
  geom_abline(slope=1, intercept=0, color = "blue") +
  xlab("Predicted Parameter") +
  ylab("True Parameter") +
  theme(legend.position = "none") +
  facet_grid(rows = vars(type), cols = vars(ntaxa)) +
  geom_text(data=myR2, size = 3, mapping = aes(x=-.5, y=3.25, label = (Rsqlab)), parse = T) +
  theme(panel.grid.minor = element_blank())
ggsave("MPDvsCKSfacet.png")
