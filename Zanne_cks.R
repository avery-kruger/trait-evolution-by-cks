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

#Read in simulated data ----
allcomm <- list.files(path="data",pattern="angiocomm", full.names = T) %>%
  map_dfr(readRDS)

#separate data into training and validation data
trainingdata <- allcomm[1:4000,]
testdata <- allcomm[4001:5000,]

Y <- trainingdata$param
X <- trainingdata$comm
valY <- testdata$param
valX <- testdata$comm

#Choose superparameters to sweep across
featuresweep <- 2^(4:9)
windowsweep <- c(5*4^(0:5),ncol(X))

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
           saveRDS(temp, paste0("data/sweep",i,".rds"))
         },
         mc.cores = numCores)

#Run this chunk if first sweep didn't finish for some reason.
#Consider skipping 29 and 30; they take a long time and are not optimal.
for(i in seq(nrow(fullsweep))[-donesweep]){
  print(i)
  temp <- kitchen_sweep(X,Y,
                        valX,valY,
                        fullsweep[i,1],fullsweep[i,2],
                        verbose = T, show.plot = F,
                        clampoutliers = T,
                        ncores = numCores)
  saveRDS(temp, paste0("data/sweep",i,".rds"))
}

#Get R2 for each combination of hyperparameters
#Alternatively, just read the file if it exists.
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

#Validation ----
##Use the best performing model to predict values for the validation data
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
mydata <- as.data.frame(mydata[[1]])
mycols <- 1:ncol(mydata)
mydata$true <- valY
mydata$mean <-sapply(seq(nrow(mydata)), function(i){
  mean(as.numeric(mydata[i,mycols]))
})
mydata$sd <-sapply(seq(nrow(mydata)), function(i){
  sd(as.numeric(mydata[i,mycols]))
})
mydata$CIlow95 <-sapply(seq(nrow(mydata)), function(i){
  quantile(as.numeric(mydata[i,mycols]), 0.025)
})
mydata$CIhigh95 <-sapply(seq(nrow(mydata)), function(i){
  quantile(as.numeric(mydata[i,mycols]), 0.975)
})

###Examine how predictions line up with true data.
ggplot(mydata,
       aes(y=true,x=mean)) +
  theme_bw() +
  geom_point(aes(alpha = 0.2)) +
  geom_errorbar(aes(xmin=CIlow95, xmax=CIhigh95))+
  geom_abline(slope=1, intercept=0, color = "blue") +
  xlab("Predicted Parameter") +
  ylab("True Parameter") +
  theme(legend.position = "none")
ggsave("ZanneValidation_bootstrap.png")

summary(lm(mydata$true ~ mydata$mean))
#R^2 is about 0.92.
sum(sign(mydata$true)==sign(mydata$mean))/nrow(mydata)
#~97.1% correct sign (will vary)

#Predict on Known Community ----
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
myalphas <- c(seq(-0.5, -0.1, 0.05), seq(-.095,-0.01,0.005),seq(-0.009, 0, 0.001),
              seq(0.001, 0.009, 0.001), seq(0.01, 0.095, 0.005), seq(0.1,0.5,0.025))
zanneMPD <- c(do.call("c", lapply(myalphas, function(x){
  mpd.query(rescale(zannetree.trim, "EB", a=x), mycomm, standardize = T)})))
names(zanneMPD) <- myalphas

###Plot MPD curve of known community ----
angiocurve <- data.frame(alpha = myalphas, MPD = zanneMPD)
ggplot(angiocurve, aes(x=alpha, y = MPD)) +
  geom_line() +
  theme_bw()
ggsave("ZanneMPDcurve.png")

###Bootstrap predictions on Zanne tree (same normal matrix)----
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
zanneboot <- zanneboot[!is.na(zanneboot)]
Zanne95mean <- mean(zanneboot)
Zanne95low <- quantile(zanneboot, 0.025)
Zanne95high <- quantile(zanneboot, 0.975)

##Main CKS Plot ----
#Plot scatter plot of predictions on validation data along with prediction of empirical data
ggplot(mydata,
       aes(y=true,x=mean)) +
  theme_bw() +
  geom_point(aes(alpha = 0.2)) +
  geom_errorbar(aes(xmin=CIlow95, xmax=CIhigh95))+
  geom_abline(slope=1, intercept=0, color = "blue") +
  xlab("Predicted Parameter") +
  ylab("True Parameter") +
  theme(legend.position = "none") +
  geom_vline(aes(xintercept = Zanne95low), col = 2, lty = 3, size = 0.5) +
  geom_vline(aes(xintercept = Zanne95mean), col = 2, lty = 2, size = 1) +
  geom_vline(aes(xintercept = Zanne95high), col = 2, lty = 3, size = 0.5) +
  xlab(expression(paste("Predicted parameter (",italic("a"),")"))) +
  ylab(expression(paste("True parameter (",italic("a"),")")))
ggsave("Zanne95boot.png")

##MPD Comparison ----
##Compare with MPD model
BigMPD <- function(BigSimList, tree, myalphas = NULL){
  if(is.null(myalphas)){
    myalphas <- myalphas <- c(seq(-0.5, -0.1, 0.05), seq(-.095,-0.01,0.005),
                              seq(-0.009, 0, 0.001), seq(0.001, 0.009, 0.001),
                              seq(0.01, 0.095, 0.005), seq(0.1,0.5,0.025))}
  MPD <- as.data.frame(cbind(do.call("cbind", lapply(myalphas, function(x){
    mpd.query(rescale(tree, "EB",x), BigSimList$comm, standardize = T)}))))
  names(MPD) <- myalphas
  MPD
}


trainingmpd <- BigMPD(trainingdata, tree = zannetree.trim)
testmpd <- BigMPD(testdata, tree = zannetree.trim)

mpd.model <- lm(Y ~ ., data=trainingmpd)
predict(mpd.model, zanneMPD.df)
mpd.val <- predict(mpd.model, testmpd)
mpd.val.model <- lm(valY ~ mpd.val)
summary(mpd.val.model)


zanneMPD.df <- as.data.frame(matrix(zanneMPD, nrow=1))
names(zanneMPD.df) <- myalphas
zanneMPD.prediction <- predict(mpd.model, zanneMPD.df)
mpd.data <- data.frame(true = valY, predict = mpd.val)

#bootstrapmodels for MPD
boot.predict <- as.data.frame(matrix(nrow = length(valY), ncol = 100))
boot.zanneMPD <- vector()
nrows <- nrow(trainingmpd)
samplestorage <- as.data.frame(matrix(nrow = 100, ncol = nrows))
set.seed(234509)
for(i in 1:100){
  mysample <- sample(seq(nrows), nrows, replace = T)
  samplestorage[i,] <- mysample
  boot.data <- trainingmpd[mysample,]
  boot.model <- lm(Y[mysample] ~ ., data=boot.data)
  boot.predict[,i] <- predict(boot.model, testmpd)
  boot.zanneMPD[i] <- predict(boot.model, zanneMPD.df)
  print(paste0(i, "/100"))
}
names(samplestorage) <- NULL

mpd.data$boot.mean <- rowMeans(boot.predict)
mpd.data$boot.95low <- sapply(seq(nrow(boot.predict)), function(i){
  quantile(as.numeric(boot.predict[i,]), 0.025)})
mpd.data$boot.95high <- sapply(seq(nrow(boot.predict)), function(i){
  quantile(as.numeric(boot.predict[i,]), 0.975)})

boot.zanneMPD.95 <- quantile(boot.zanneMPD, c(0.025,0.975))

####Plot MPD predictions ----
ggplot(mpd.data,
       aes(y=true,x=predict)) +
  theme_bw() +
  geom_point(aes(alpha = 0.2)) +
  geom_errorbar(aes(xmin=boot.95low, xmax=boot.95high))+
  geom_abline(slope=1, intercept=0, color = "blue") +
  xlab("Predicted Parameter") +
  ylab("True Parameter") +
  theme(legend.position = "none") +
  geom_vline(aes(xintercept = boot.zanneMPD.95[1]), col = 2, lty = 3, size = 0.5) +
  geom_vline(aes(xintercept = zanneMPD.prediction), col = 2, lty = 2, size = 1) +
  geom_vline(aes(xintercept = boot.zanneMPD.95[2]), col = 2, lty = 3, size = 0.5) +
  xlab(expression(paste("Predicted parameter (",italic("a"),")"))) +
  ylab(expression(paste("True parameter (",italic("a"),")")))
ggsave("MPDboot_95CI.png")




###Rescaled Tree with Color----
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
