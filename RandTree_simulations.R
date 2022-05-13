library(geiger)
library(picante)
library(phytools)
library(PhyloMeasures)
library(ape)
library(reshape2)
library(parallel)
source("sim_MPD_functions.R")

MPDapply <- function(community, tree=mytree, myalphas = NULL){
  myalphas <- c(seq(-3, -0.1, 0.1), seq(-.095,-0.01,0.005),seq(-0.009, 0, 0.001),
                seq(0.001, 0.009, 0.001), seq(0.01, 0.095, 0.005), seq(0.1,3,0.1))
  MPD <- do.call("cbind", lapply(myalphas, function(x){
    mpd.query(rescale(tree, "EB", a=x), community, standardize = T)}))
  colnames(MPD) <- myalphas
  MPD
}

community_sim <- function(tree,
                          reps=100,
                          comm_size = 30,
                          paramFUN = `rnorm`,
                          param_list = list(n=1,mean=0,sd=0.5),
                          param_bounds = NULL,
                          optimum_species = NULL,
                          verbose = TRUE){
  if(length(tree$tip.label) < comm_size){stop('Tree larger than community size')}
  mytree <- tree
  n.sp <- length(mytree$tip.label)

  Communities <- list()
  Communities$comm <- as.data.frame(matrix(ncol = length(mytree$tip.label)))
  colnames(Communities$comm) <- mytree$tip.label
  Communities$filter <- NULL #for storing filter strength (1 = no sp)
  Communities$param <- NULL #for storing tree transformation parameter

  abundances <- 1
  J <- comm_size#filter strength
  blankcomm <- rep(0,n.sp)
  time <- Sys.time()
  for(i in 1:reps){
    param <- do.call(paramFUN, param_list)
    if(!is.null(param_bounds)){
      param <- min(max(param, param_bounds[1]),param_bounds[2])
    }
    newtree <- rescale(mytree, model = "EB", a=param)

    bm<-fastBM(newtree, nsim=10)

    if(is.null(optimum_species)){
      filter.optimum <- bm[sample(1:n.sp,1),]
    } else{
      filter.optimum <- bm[sample(optimum_species,1),]
    }
    head(filter.optimum)
    mydist <- traitdist(bm, filter.optimum)

    mycomm <- closestN(N = comm_size,mydist,1)


    Communities$comm[i,] <- mycomm
    Communities$param[i] <- param
    if(Sys.time() - time > 0.5 & verbose){
      print(paste0("Simulation ",i,"/",reps))
      time <- Sys.time()
    }
  }
  Communities
}


set.seed(371015)
for(ntaxa in c(48, 2^(6:10))){
  mytree <- drop.extinct(sim.bdtree(.8,0.2,"taxa",ntaxa))
  commsize <- 30
  filename <- paste0(ntaxa,"treesim")
  for(i in 1:80){
    newsim <- community_sim(mytree,
                                   reps=25,
                                   comm_size = commsize,
                                   param_list = list(n=1,mean=0,
                                                     sd=5/max(node.depth.edgelength(mytree)))

    )
    saveRDS(newsim, file = paste0("data/",filename,"_comm",i,".rds"))
    print(paste0(ntaxa,"tree: ",i,"/",80, " Community Simulated"))
    mpd <- list()
    mpd$mpd <- MPDapply(newsim$comm, tree=mytree)
    mpd$param <- newsim$param
    saveRDS(mpd, file = paste0("data/",filename,"_mpd",i,".rds"))
    print(paste0(ntaxa,"tree: ",i,"/",80, " MPD Calculated"))
  }
}


