#gets Euclidean trait distances between all species and some optimum.
traitdist <- function(traits, optimum){
  sqrt(rowSums(sweep(traits, 2, optimum)^2))
}

#gets the closest n species given their distances from an optimum
closestN <- function(N, traitdist, prob = 1){ #prob is just for working with coalesc
  orderdist <- sort(unique(traitdist)) #[1:N]
  chosen <- orderdist[orderdist < orderdist[N]]
  remainder <- N - length(chosen)
  newlist <- traitdist
  newlist[which(traitdist %in% chosen)] <- prob
  #when there is a tie (because of extremely similar values, low or high)
  #it randomly selects from the tied species to fill the remainder.
  if(remainder > 1){
    newlist[sample(which(traitdist == orderdist[N]),remainder)] <- prob
  } else{
    newlist[which(traitdist == orderdist[N])] <- prob}
  newlist[newlist != prob] <- 0
  newlist
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
    #param <- runif(1, paramrange[1], paramrange[2]) #too few possible sp sometimes when param gets to around >.9
    newtree <- rescale(mytree, model = "EB", a=param)

    bm<-fastBM(newtree, nsim=10)

    if(is.null(optimum_species)){
      filter.optimum <- bm[sample(1:n.sp,1),]
    } else{
      filter.optimum <- bm[sample(optimum_species,1),]
    }
    head(filter.optimum)
    mydist <- traitdist(bm, filter.optimum)
    print(head(mydist))
    print(param)
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

MPDapply <- function(communities, tree=mytree, mydeltas = NULL){
  if(is.null(mydeltas)){mydeltas <- c(seq(0.05,0.95,.05),seq(1,40,0.5))}
  MPD <- as.data.frame(cbind(do.call("cbind", lapply(mydeltas, function(x){
    mpd.query(rescale(tree, "delta",x), communities, standardize = T)}))))
  names(MPD) <- mydeltas
  MPD
}
