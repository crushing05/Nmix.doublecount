#' sim_data
#'
#' Simulate counts and capture-recapture data
#' 
#' @param S Number of sites
#' @param K Number of visits
#' @param lambda Expected abundance 
#' @param phi Proportion of individuals that are marked
#' @param theta Visitation rate
#' @export

sim_data <- function(S, K, lambda, phi, theta) {

  ## Abundance
  N <- rpois(S, lambda)                    # Site-level abundance
  Ntot <- sum(N)                           # Total abundance
  site <- rep(1:S, N)                      # Site associated with each individual.
  
  
  ## Marked status
  marked <- rbinom(Ntot, 1, phi)           # Is individual marked (Y/N)
  TotNm <- sum(marked)                     # Total number of marked individuals
  markedSite <- site[which(marked==1)]     # Site associated with each marked individual
  unmarkedSite<-site[which(marked==0)]     # Site associated with each unmarked individual
   
  ## Marked data
  Nm <- rep(0, S)                      # Placehoder for marked individuals per site
  Nu <- rep(0, S)                      # Placehoder for unmarked individuals per site
  Nm[as.numeric(names(table(markedSite)))] <- table(markedSite)       # Number of marked individuals per site
  Nu[as.numeric(names(table(unmarkedSite)))] <- table(unmarkedSite)   # Number of unmarked individuals per site
  
  ## True number of detections of each individual
  y.true <- matrix(rpois(Ntot*K, theta), nrow = Ntot, ncol = K) 
  
  ## Observed number of detections of marked individuals
  y <- y.true[which(marked == 1),]     
  
  ## Counts of unmarked individuals
  y.unmarked <- y.true[which(marked == 0),]     # Number of detections of unmarked individuals
  n <- n.all <- matrix(NA, nrow = S, ncol = K)  # Matrices to store counts of marked and all individuals 
  
  for(i in 1:S){
    n[i,] <- apply(y.unmarked[unmarkedSite == i,], 2, sum)
    n.all[i,] <- apply(y.true[site == i,], 2, sum)
  }
  
  dat <- list(K = K, S = S, N = N, lambda = lambda, phi = phi, theta = theta, 
              Ntot = Ntot, TotNm = TotNm, Nm = Nm, 
              Nu = Nu, y.true = y.true, y = y, n = n, n.all = n.all)
  return(dat)
}
