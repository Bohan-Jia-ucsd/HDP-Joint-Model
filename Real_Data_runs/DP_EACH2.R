SEED.START = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
library(tidyverse)
library(R2jags)
load(paste0("dat_CD4_combine.rdata"))
load(paste0("dat_surv_combine.rdata"))

dat_surv_combine<-dat_surv_combine%>%filter(pop=="320")
dat_CD4_combine<-dat_CD4_combine%>%filter(pop=="320")
n=nrow(dat_surv_combine)
dat_surv_combine<-dat_surv_combine%>%mutate(id=1:n)
remained_id<-dat_CD4_combine%>%group_by(id)%>%summarise(n=n())
id_new<-c()
for (i in 1:n) {
  id_new<-c(id_new,rep(i,remained_id$n[i]))
}
dat_CD4_combine$id<-id_new

# model
cat("model{
  
  for (k in 1:K)
  {
    ###survival
    m[k] ~ dnorm(0,0.1) # parameters for y
    
    ###lmm
    intercept[k] ~ dnorm(0,0.1) # fixed intercept
    tau_b[k] ~ dgamma(1,0.1)
    sigma_b[k] <- sqrt(1/tau_b[k]) # sd for random intercept
  }
  
    tau ~ dgamma(1,0.1)
    sigma <- sqrt(1/tau)
    a ~ dnorm(0,0.1)
    base_haz~ dgamma(1,0.1)
  
  for (i in 1:n){b[i] ~ dnorm(0,tau_b[zeta[i]])}
  
  for (i in 1:N){
    mu[i]<-a*t[i]+b[id[i]]+intercept[zeta[id[i]]]
    y[i] ~ dnorm(mu[i],tau)
  }
  
  for (i in 1:n){
    y_estimated[i]<-Surv_t[i]*a+b[id_Surv[i]]+intercept[zeta[id_Surv[i]]]
    
    elinpred[i] <- exp(m[zeta[id_Surv[i]]]*y_estimated[i])
    
    logHaz[i] <-log(base_haz*(Surv_t[i]^(base_haz-1))*elinpred[i])
    
    ### integral of the survival function
    for(l in 1:N_grid){
      f[i,l]<-base_haz*(((Surv_t_grid_point[i,l]+Surv_t_grid_point[i,l+1])/2)^(base_haz-1))*
        exp(m[zeta[id_Surv[i]]]*
              (b[id_Surv[i]]+intercept[zeta[id_Surv[i]]]))*
        exp(m[zeta[id_Surv[i]]]*a*((Surv_t_grid_point[i,l]+Surv_t_grid_point[i,l+1])/2))*
        (Surv_t_grid_point[i,l+1]-Surv_t_grid_point[i,l])}
    
    logSurv[i] <- -sum(f[i,1:N_grid])
    
    phi[i]<- 100 - delta[i]*logHaz[i]-logSurv[i]
    zeros[i]~dpois(phi[i])
  }
  
 for (i in 1:n){zeta[i] ~ dcat(pi[])}
  
 # Stick breaking
 rho ~ dgamma(1,1)
 for (k in 1:(K-1)) { V[k] ~ dbeta(1,rho)T(0.0001,0.9999)}
 V[K] <- 1
 pi[1] <- V[1]
 
 for (k in 2:K) {
 pi[k] <- V[k] * (1-V[k-1]) * pi[k-1]/V[k-1]
 }
  
  
}
",file="DP_each2_joint.txt")

K=10
N_grid=100

Surv_t_grid_point<-matrix(0,nrow =nrow(dat_surv_combine), ncol = N_grid+1)
for (i in 1:nrow(dat_surv_combine)) {
  Surv_t_grid_point[i,1:(N_grid+1)]<-seq(from=0, to=dat_surv_combine$ART_init[i], 
                                         by=dat_surv_combine$ART_init[i]/N_grid)
}
n=nrow(dat_surv_combine)
data_test<-list(K=K,
                n=nrow(dat_surv_combine),
                #lmm
                id=dat_CD4_combine$id,
                N=nrow(dat_CD4_combine),
                t=dat_CD4_combine$visit,
                y=dat_CD4_combine$cd4count,
                #surv
                Surv_t=dat_surv_combine$ART_init,
                id_Surv=dat_surv_combine$id,
                zeros=rep(0,nrow(dat_surv_combine)),
                delta=dat_surv_combine$status,
                #integral
                N_grid=N_grid,
                Surv_t_grid_point=Surv_t_grid_point)
set.seed(SEED.START^2+200*SEED.START)
jags.m <- jags( model.file= "DP_each2_joint.txt", 
                data=data_test, 
                inits=list(list(rho=rgamma(1,1,1),
                                beta=rnorm(K,0,1),
                                m=rnorm(K,0,1),base_haz=rgamma(1,1,1),
                                tau_b=rgamma(K,1,0.1),a=rnorm(1,0,1),
                                tau=rgamma(1,1,0.1),intercept=rnorm(K,0,1))),
                n.chains=1, 
                n.iter = 20000,
                n.burnin =17000,
                n.thin=1,
                parameters.to.save=c("base_haz","beta","m","pi",
                                     "sigma_b","sigma","intercept","a","zeta")
)

save(jags.m, file = paste0("DP_Each2_594320.rdata"))



