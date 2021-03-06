sink("inst/jags/Nmix_binom.jags")
cat("
    model {
    ##################
    ### Priors
    ##################

      lambda ~ dunif(0, 100)         # Abundance rate
      p ~ dunif(0, 1)                # Detection probability
    

    ##################
    ### Likelihoods
    ##################

      for (s in 1:S) {
        ## Total abundance
        N[s] ~ dpois(lambda)
      
      
        ## Counts
        for (k in 1:K) {
          n[s,k] ~ dbinom(p, N[s])           
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