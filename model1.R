##########################################################################
##########################################################################
##This script runs a hierarchical site occupancy-detection model in JAGS##
#####It uses a 4D array  [visit, site, sample,qPCR] and covarites for#####
##each time period post fire, burnt status of the site and water volume###
########################for each sample###################################
##########################################################################
##########################################################################


# packages ----------------------------------------------------------------

library(jagsUI)



# data --------------------------------------------------------------------

load("data.RData")

#model1_sp_data= 4d array [visit, site, sample,qPCR]
#model1_I_after1=1 year post-fire
#model1_I_after2=2 years post-fire
#model1_I_burnt=Burnt status of watershed
#model1_volume=volume of water filtered per sample


#model -----------------------------------------------------------


str(model1<-list(model1_sp_data = model1_sp_data,nvisit=dim(model1_sp_data)[1], nsite = dim(model1_sp_data)[2], 
                 nrep = dim(model1_sp_data)[3], 
                 npcr = dim(model1_sp_data)[4],
                 I_after1=model1_I_after1,
                 I_after2=model1_I_after2,
                 I_burnt=model1_I_burnt,
                 volume=model1_volume))

sink("model1.txt")
cat("
    model {
    
    # Priors
    
    # - Occupancy
    for(j in 1:nsite){
    int.psi[j] ~ dunif(0,1)
    }
    # - eDNA availability
    
    int.theta ~ dunif(0,1) # Intercepts availability probability
    
    #detection
    int.p ~ dunif(0,1)
  
    
    # Hyperpriors
    # - occupancy
    betaBA1 ~ dnorm(0, 0.1)
    betaBA1burnt ~ dnorm(0, 0.1)
    betaBA2 ~ dnorm(0, 0.1)
    betaBA2burnt ~ dnorm(0, 0.1)
    betaBurnt ~ dnorm(0, 0.1)
    beta_volume~ dnorm(0, 0.1)
    
    
    # Ecological process 
    for(i in 1:nvisit){
    for(j in 1:nsite){
    logit(pocc[i, j]) <-logit(int.psi[j])  + betaBA1 * I_after1[i, j]  + betaBA2 * I_after2[i, j] +betaBurnt* I_burnt[i, j]+
    betaBA1burnt * I_after1[i, j] * I_burnt[i, j] + betaBA2burnt * I_after2[i, j] * I_burnt[i, j]
    z[i, j] ~ dbern(pocc[i, j])
    
    
    for (k in 1:nrep){
    # Occurrence in sample j
    a[i,j,k] ~ dbern(mu.a[i,j,k])
    mu.a[i,j,k] <- z[i, j] * theta[i,j,k]
    cloglog(theta[i,j,k]) <- cloglog(int.theta)+beta_volume*volume[i,j,k]
    
    
    
    # Observation process
    
    for(l in 1:npcr){
    
    logit(p[i, j, k, l])<- logit(int.p)
    mu.p[i, j, k, l] <-  a[i,j,k] *p[i, j,k,l]
    
    model1_sp_data[i, j, k, l] ~ dbern(mu.p[i, j, k,l])
    
    }
    }
    }
    }
    
    }
    ",fill = TRUE)
sink()

zst<-matrix(data=1,nrow=3,ncol= dim(model1_sp_data)[2])
ast<-array(1, c(3, dim(model1_sp_data)[2], 2))

ast <- apply(model1_sp_data, c(1,2,3), max)   # inits for availability (a)
ast[is.na(ast)] <- 1

inits <- function() list(z = zst, a = ast) 
# Parameters monitored 
params1 <- c("int.psi","int.theta", "int.p",
             "betaBA1","betaBA2","betaBA1burnt",
             "betaBA2burnt","betaBurnt", "beta_volume","z",
             "pocc") 


# MCMC settings
ni <- 30000    ;   nt <- 10   ;   nb <- 10000   ;   nc <- 3


out_model1<- jags(model1, inits, params1, "model1.txt", n.chains = nc, n.thin = nt, n.iter = ni,n.burnin = nb)

print(out_model1, 3)

