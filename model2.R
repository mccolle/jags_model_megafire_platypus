
##########################################################################
#########################################################################
#this script formats the presence/absence data (sp_data) and covariates##
#into the correct format for jags model##################################
#



# packages ----------------------------------------------------------------

library(jagsUI)



# data --------------------------------------------------------------------



load("model2_sp_data.RData") #4d array [visit, site, sample,qPCR]
load("model2_I_after1.RData") #1 year post-fire
load("model2_I_after2.RData") #2 years post-fire
load("model2_h_sev.RData") #high severity fire
load("model2_l_sev.RData") #low severity fire
load("model2_rain.RData") #rainfall
load("model2_volume.RData") #volume of water filtered per sample




# model -------------------------------------------------------------------


str(model2<-list(sp_data = model2_sp_data,nvisit=dim(model2_sp_data)[1], nsite = dim(model2_sp_data)[2], 
                 nrep = dim(model2_sp_data)[3], 
                 npcr = dim(model2_sp_data)[4],
                 I_after1=model2_I_after1,
                 I_after2=model2_I_after2,
                 h_sev=rbind(model2_h_sev,model2_h_sev,model2_h_sev),
                 l_sev=rbind(model2_l_sev,model2_l_sev,model2_l_sev),
                 rain=rbind(model2_rain,model2_rain,model2_rain),
                 volume=model2_volume
))


sink("model2.txt")
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
    
    
    
    # - Before-after effect
    betaBA1 ~ dnorm(0, 0.1)
    betaBA2 ~ dnorm(0, 0.1)
    
    
    #covs
    beta_lsev~ dnorm(0, 0.1)
    beta_hsev~ dnorm(0, 0.1)
    beta_rain ~ dnorm(0,0.1)
    beta_volume~ dnorm(0,0.1)
    
    #two way interactions
    
    beta_i_rain_hsev~ dnorm(0, 0.1)
    beta_i_rain_lsev ~ dnorm(0, 0.1)
    
    beta_lsev_BA1~ dnorm(0, 0.1)
    beta_lsev_BA2~ dnorm(0, 0.1)
    beta_hsev_BA1~ dnorm(0, 0.1)
    beta_hsev_BA2~ dnorm(0, 0.1)
    beta_rain_BA1~ dnorm(0, 0.1)
    beta_rain_BA2~ dnorm(0, 0.1)

    
    # Ecological process 
    for(i in 1:nvisit){
    for(j in 1:nsite){
    logit(pocc[i, j]) <-
    logit(int.psi[j]) + betaBA1 * I_after1[i, j]  +
    betaBA2 * I_after2[i,j] + 
    beta_hsev*h_sev[i,j]+
    beta_lsev*l_sev[i,j]+
    beta_rain*rain[i,j]+
    
    
    beta_i_rain_hsev*rain[i,j]*h_sev[i,j]+
    beta_i_rain_lsev*rain[i,j]*l_sev[i,j]+
    
    beta_lsev_BA1*l_sev[i,j]* I_after1[i, j]+
    beta_lsev_BA2*l_sev[i,j]* I_after2[i, j]+
    beta_hsev_BA1*h_sev[i,j]* I_after1[i, j]+
    beta_hsev_BA2*h_sev[i,j]* I_after2[i, j]+
    
    beta_rain_BA1*rain[i,j]* I_after1[i, j]+
    beta_rain_BA2*rain[i,j]* I_after2[i, j]
    
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
    sp_data[i, j, k, l] ~ dbern(mu.p[i, j, k,l])
    }
    }
    }
    }
    }
    ",fill = TRUE)
sink()

zst<-matrix(data=1,nrow=3,ncol=nsite)
ast<-array(1, c(3,nsite, 2))

ast <- apply(sp_data, c(1,2,3), max)   # inits for availability (a)
ast[is.na(ast)] <- 1

inits <- function() list(z = zst, a = ast)

# Parameters monitored 
params1 <- c("int.psi","int.theta", "int.p",
             "betaBA1","beta_hsev",
             "betaBA2",
             "beta_lsev","beta_rain","beta_volume",
             "beta_i_rain_hsev","beta_i_rain_lsev",
             "beta_lsev_BA1","beta_lsev_BA2","beta_hsev_BA1","beta_hsev_BA2",
             "beta_rain_BA1","beta_rain_BA2") 


# MCMC settings
ni <- 30000    ;   nt <- 10   ;   nb <- 10000   ;   nc <- 3

out_model2<- jags(model2, inits, params1, "model2.txt", n.chains = nc, n.thin = nt, n.iter = ni,n.burnin = nb)
print(out_model2, 3)
