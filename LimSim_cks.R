library(geiger)
library(picante)
library(phytools)
library(PhyloMeasures)
library(ecolottery)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(visreg)
library(ggtree)
library(kitchen)
source("sim_MPD_functions.R")

allsims <- list.files(path="data",pattern="limsim_sim", full.names = T) %>%
  map_dfr(readRDS)

allcomm <- list.files(path="data",pattern="limsim_comm", full.names = T) %>%
  map_dfr(readRDS)

trainingdata <- allcomm[1:1200,]
testdata <- allcomm[1201:2000,]

Y <- trainingdata$param
X <- trainingdata$comm #[,-which(names(trainingdata) == "param")]
valY <- testdata$param
valX <- testdata$comm #[,-which(names(trainingdata) == "param")]

#Choose superparameters to sweep across
featuresweep <- 2^(4:9)
windowsweep <- c(5, 10, 20, 40, 80, 128) #c(5*4^(0:5),ncol(X)) #c(5,10,20,40,98)

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
ggsave("LimSim_128SpTree.png")

summary(lm(mydata$true ~ mydata$predicted))
#R^2 is about 0.85.
sum(sign(mydata$true)==sign(mydata$predicted))/nrow(mydata)
#~97.7% correct sign (will vary)

#Compare with MPD
trainingmpd <- allsims[1:1200,]$mpd
testmpd <- allsims[1201:2000,]$mpd
trainingYmpd <- allsims[1:1200,]$param
testYmpd <- allsims[1201:2000,]$param

mpd.model <- lm(trainingYmpd ~ ., data=as.data.frame(trainingmpd))
mpd.val <- predict(mpd.model, as.data.frame(testmpd))
mpd.val.model <- lm(testYmpd ~ mpd.val)
summary(mpd.val.model)
par(xpd=FALSE)
plot(mpd.val, testYmpd)
abline(0, 1)
MPDdata <- data.frame(true = testYmpd, predicted = mpd.val)

ggplot(MPDdata,
       aes(y=true,x=predicted)) +
  theme_bw() +
  geom_point(aes(alpha = 0.2)) +
  geom_abline(slope=1, intercept=0, color = "blue") +
  xlab("Predicted Parameter") +
  ylab("True Parameter") +
  theme(legend.position = "none")
ggsave("LimSim_128SpTree_mpd.png")

