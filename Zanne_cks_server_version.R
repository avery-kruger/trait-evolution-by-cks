library(geiger)
library(picante)
library(phytools)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(ggtree)
library(devtools)
library(parallel)
library("kitchen")
#source("sim_MPD_functions.R")

allcomm <- readRDS("data/allcomm.rds")

#separate data into training and validation data
trainingdata <- allcomm[1:4000,]
testdata <- allcomm[4001:5000,]

Y <- trainingdata$param
X <- trainingdata$comm #[,-which(names(trainingdata) == "param")]
valY <- testdata$param
valX <- testdata$comm #[,-which(names(trainingdata) == "param")]

#Choose superparameters to sweep across
featuresweep <- 2^(4:9)
windowsweep <- c(4^(2:6),ncol(X))#c(5*4^(0:5),ncol(X)) #c(5,10,20,40,98)
fullsweep <- expand.grid(featuresweep,windowsweep)
#For running only part, or incase of crashes, saves progress.
donesweep <- as.numeric(
  gsub("sweep([0-9]+).rds","\\1",list.files(path="data",pattern="sweep")))

#Conduct preliminary sweep
numCores <- detectCores()
mclapply(seq(nrow(fullsweep)),
         function(i){
           temp <- kitchen_sweep(X,Y,
                                 valX,valY,
                                 fullsweep[i,1],fullsweep[i,2],
                                 verbose = T, show.plot = F)
           saveRDS(temp, paste0("sweep",i,".rds"))
         },
         mc.cores = numCores)
#saveRDS(mysweep,"mysweep.rds")

#Consider skipping 29 and 30; they take a long time and are not optimal.
for(i in seq(nrow(fullsweep))[-donesweep]){
  print(i)
  temp <- kitchen_sweep(X,Y,
                        valX,valY,
                        fullsweep[i,1],fullsweep[i,2],
                        verbose = T, show.plot = F,
                        clampoutliers = T,
                        ncores = numCores)
  saveRDS(temp, paste0("sweep",i,".rds"))
}

#Get R2 for each combination of hyperparameters
#Alternatively, just read the file if you have it.
if(length(list.files(path="data","fullsweepR2.rds")) == 0){
  allsweep_names <- paste0("data/",list.files(path = "data", pattern = "sweep"))
  allsweep_index <- as.numeric(
    gsub("data/sweep([0-9]+).rds","\\1",allsweep_names))
  allsweep <- list()
  fullsweep$r2 <- NA
  for(i in seq(length(allsweep_names))){
    allsweep[[i]] <- readRDS(paste0(allsweep_names[[i]]))
    fullsweep$r2[allsweep_index[i]] <- as.numeric(allsweep[[i]])
  }
  saveRDS(fullsweep, "data/fullsweepR2.rds")
}

fullsweep <- readRDS("fullsweepR2.rds")

#bestR2
best.index <- which(fullsweep$r2 == max(fullsweep$r2,na.rm=T))
fullsweep[best.index, 1] #best feature
fullsweep[best.index, 2] #best window

###Use the best performing model to predict values for the validation data
if("validation_bootstrap100.rds" %in% list.files()){
  mydata <- readRDS("data/validation_bootstrap100.rds")
} else{mydata <- data.frame(true = valY)}
if(ncol(mydata) < 100){
  for(i in ncol(mydata):100){
    bestmodelpredictions <- kitchen_prediction(X,
                                               Y,
                                               valX,
                                               fullsweep[best.index, 1],
                                               fullsweep[best.index, 2],
                                               verbose = T,
                                               bootstrap = 100,
                                               write_progress = "data/validation_bootstrap100.rds",
                                               ncores = numCores)

    mydata[,paste0("predicted",i)] <- bestmodelpredictions[[1]][[1]]
    saveRDS(mydata, "data/validation_bootstrap100.rds")
    print(paste0(i,"/100 done"))
  }
}

mydata <- readRDS("data/validation_bootstrap100.rds")

###Examine how predictions line up with true data.
mydata$true <- valY
ggplot(mydata,aes(x=mean, y=true)) +
  theme_bw() +
  geom_point() +
  xlab("Predicted Parameter") +
  ylab("True Parameter") #+
  #geom_vline(xintercept = meanpred) +
  #geom_abline(0, 1)
ggsave("ZanneValidation.png")

summary(lm(mydata$true ~ mydata$mean))
#R^2 is about 0.85.
sum(sign(mydata$true)==sign(mydata$mean))/nrow(mydata)
#~97.7% correct sign (will vary)

#is this necessary? maybe remove
#Compare with simple lm
# mylm <- lm(Y ~ ., data=X)
# summary(mylm)
# lmpred <- predict(mylm, valX)
# summary(lm(valY ~ lmpred))
# plot(valY ~ lmpred)

###Test w/ real Zanne data
zannetree <- read.tree("Zanne.angiosperm.tre")
zannefreeze <- read.csv("MinimumFreezingExposure.csv")
zannetree$tip.label <- trimws(tolower(zannetree$tip.label),"both")
zannefreeze$species <- trimws(gsub(" ","_",tolower(zannefreeze$species)),"both")
nofreezenames <- which(is.na(match(zannetree$tip.label,zannefreeze$species)))
myfreeze <- zannefreeze[which(!is.na(match(zannefreeze$species,zannetree$tip.label))),]
zannetree.trim <- drop.tip(zannetree,nofreezenames)
commsize <- sum(myfreeze$Freeze.tmin.lo == "FreezingExposed")

mycomm <- 1*(myfreeze$Freeze.tmin.lo == "FreezingExposed")
mycomm <- t(as.data.frame(mycomm, row.names = myfreeze$species))

mynorm <- make_norms(fullsweep[best.index, 1],fullsweep[best.index, 2])[[1]][[1]]


zanneboot <- kitchen_prediction(X,
                                Y,
                                mycomm,
                                fullsweep[best.index, 1],
                                fullsweep[best.index, 2],
                                verbose = T,
                                simplify = F,
                                bootstrap = 100,
                                write_progress = "zanneboot.rds",
                                seed = 234587,
                                ncores = numCores)


zanneboot <- unlist(readRDS("data/zanneboot.rds"))
zanneboot <- c(unlist(zanneboot),unlist(zanneboot2),unlist(zanneboot3))
zanneboot <- zanneboot[!is.na(zanneboot)]
meanpred <- mean(unlist(zanneboot), na.rm=T)
zanneCF <- sort(zanneboot)[c(floor(.025*length(zanneboot)),floor(.975*length(zanneboot)))]

mypred <- fullsweep[,1:2]
mypred$pred <- NA
for(i in 1:23){
  mypred$pred[i] <- unlist(readRDS(paste0("data/zannepred",i,".rds")))
}
mean(mypred$pred, na.rm=T)


mydeltas <- c(seq(0.05,0.95,.05),seq(1,40,0.5))
myalphas <- c(seq(-0.5, -0.1, 0.05), seq(-.095,-0.01,0.005),seq(-0.009, 0, 0.001),
              seq(0.001, 0.009, 0.001), seq(0.01, 0.095, 0.005), seq(0.1,0.5,0.025))
zanneMPD <- c(do.call("c", lapply(myalphas, function(x){
  mpd.query(rescale(zannetree.trim, "EB", a=x), mycomm, standardize = T)})))
names(zanneMPD) <- myalphas

angiocurve <- data.frame(alpha = myalphas, MPD = zanneMPD)
ggplot(angiocurve, aes(x=alpha, y = MPD)) +
  geom_line() +
  theme_bw()
ggsave("ZanneMPDcurve.png")

###Bootstrap predictions with 5,000 unique normal matrices using the same
###superparameters.
Zanne_prediction <- kitchen_prediction(X,Y,
                                       zanneMPD,
                                       featuresweep[best.index[1]],
                                       windowsweep[best.index[2]],
                                       reps=5000, simplify = T)

hist(Zanne_prediction,breaks=40)
meanpred <- mean(Zanne_prediction)
ZannePrediction.95CI <- sort(Zanne_prediction)[c(.025*5000,.975*5000)]
Zanne95low <- ZannePrediction.95CI[1]
Zanne95high <- ZannePrediction.95CI[2]

#Plot one example of truth vs prediction with 95% confidence interval from
#known data overlay.
ggplot(mydata,
       aes(y=true,x=predicted)) +
  theme_bw() +
  geom_point() +
  geom_vline(aes(xintercept = zanneboot[2]), col = 2, lty = 3, size = .5) +
  geom_vline(aes(xintercept = mean(zanneboot)), col = 2, lty = 2, size = 1) +
  geom_vline(aes(xintercept = zanneboot[25]), col = 2, lty = 3, size = .5) +
  xlab(expression(paste("Predicted parameter (",italic("a"),")"))) +
  ylab(expression(paste("True parameter (",italic("a"),")")))
ggsave("Zanne95CI.png")


##Compare with MPD model
BigMPD <- function(BigSimList, tree, mydeltas = NULL){
  if(is.null(mydeltas)){mydeltas <- c(seq(0.05,0.95,.05),seq(1,40,0.5))}
  MPD <- as.data.frame(cbind(do.call("cbind", lapply(mydeltas, function(x){
    mpd.query(rescale(tree, "delta",x), BigSimList$comm, standardize = T)}))))
  names(MPD) <- mydeltas
  MPD
}

trainingmpd <- BigMPD(trainingdata, tree = zannetree.trim)
testmpd <- BigMPD(testdata, tree = zannetree.trim)
mpd.model <- lm(Y ~ ., data=trainingmpd)
mpd.val <- predict(mpd.model, testmpd)
mpd.val <- clamp(mpd.val, min(valY))
mpd.val.model <- lm(valY ~ mpd.val)
summary(mpd.val.model)
par(xpd=FALSE)
plot(mpd.val, valY)
abline(0, 1)
zanneMPD.df <- as.data.frame(matrix(zanneMPD, nrow=1))
names(zanneMPD.df) <- mydeltas
zanneMPD.prediction <- predict(mpd.model, zanneMPD.df)

###Plot tree, colored and with scale bars according to transformation
treescale <- function(breakpoints, labels = NULL, y=1, tickheight=20,alternate = F, ...){
  #print(length(breakpoints))
  for(i in 1:(length(breakpoints)-1)){
    myheight <- tickheight
    if(alternate && i%%2 == 0){myheight <- -tickheight}
    segments(breakpoints[i],y,breakpoints[i],y+myheight)
    segments(breakpoints[i],y,breakpoints[i+1],y)
    segments(breakpoints[i+1],y,breakpoints[i+1],y+myheight)
    if(!is.null(labels)){text(breakpoints[i],y+4*myheight,labels[i],...)}
  }
  if(alternate && length(breakpoints)%%2 == 0){myheight <- -tickheight}
  if(!is.null(labels)){text(breakpoints[length(breakpoints)],y+4*myheight,labels[length(breakpoints)],...)}
}

zp <- meanpred
nticks <- 5
maxage <- max(node.depth.edgelength(zannetree.trim))
ticks <- c(seq(0,maxage,by=100),maxage)[-5]
scaleticks <- c(ticks[1],log(ticks[2:(length(ticks)-1)])/zp,ticks[length(ticks)])
scaleticks <- max(ticks)-exp(ticks*zp)/max(exp(ticks*zp))*max(ticks)
mylabels <- round(rev(ticks),0)
ticks2 <- c(seq(0,maxage,by=100),maxage)[-5]

#r[t] = 1 * exp(a * t)
#integral of rate from time t0 to time t1 = exp(a * t) - e^0 = exp(at1) -1

myratio <- rescale(zannetree.trim, "EB", zp)$edge.length/zannetree.trim$edge.length
myratio <- myratio/(max(myratio))
treecolor <- colorRampPalette(c("black","red"))
mycolors <- treecolor(51)[round(50*myratio + 1)]

op <- par
png("ZanneTreeTransformScale.png", width = 1200, height = 900, res = 200, units = "px")
par(xpd=T, mar=par()$mar+c(1,0,0,0))
plot(zannetree.trim, show.tip.label = FALSE, edge.color = mycolors)
treescale(rev(round(scaleticks,0)),y=-1000,tickheight = 150,labels = (round(mylabels,0)), cex=0.5, alternate = T)
treescale(rev(round(abs(ticks2-max(ticks2)))),y=-2000,tickheight = -150,labels = rev(round(ticks2,0)), cex=0.5)
text(7, -455, "MA (Transformed)", cex = 0.5, adj = c(0,NA))
text(7, -2590, "MA", cex = 0.5, adj = c(0,NA))
dev.off()
