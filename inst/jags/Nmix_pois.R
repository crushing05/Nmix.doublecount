sink("inst/jags/Nmix_pois.jags")
cat("
    model {
    ##################
    ### Priors
    ##################

      lambda ~ dunif(0, 100)         # Abundance rate
      theta ~ dunif(0, 100)          # Detection rate
    

    ##################
    ### Likelihoods
    ##################

      for (s in 1:S) {
        ## Total abundance
        N[s] ~ dpois(lambda)
      
      
        ## Counts
        for (k in 1:K) {
          n[s,k] ~ dpois(theta * N[s])         
        }#k
      }#s
    

    ###########################
    ### Derived parameters
    ###########################
    
    ## Total abundance

      Ntot <- sum(N[1:S])

    }#model
    ",fill = TRUE)
sink()