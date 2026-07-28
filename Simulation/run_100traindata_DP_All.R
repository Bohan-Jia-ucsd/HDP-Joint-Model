
SEED.START = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
library(tidyverse)
library(R2jags)
load( paste0("dat_train",SEED.START,".rdata"))
dat_lmm<-dat%>%select(-event_time,-status)
dat_surv<-dat%>%group_by(id)%>%slice(1)%>%ungroup()%>%
  select(-y,-t)
n_aid<-nrow(dat_surv[dat_surv$J==1,])

# model
cat("model{
  
  for (k in 1:K)
  {
    ###survival
    m[k] ~ dnorm(0,0.1) # parameters for y
    beta[k] ~ dnorm(0,0.1) # parameters for baseline covariates
    
    ###lmm
    intercept[k] ~ dnorm(0,0.1) # fixed intercept
    tau_b[k] ~ dgamma(1,0.1)
    sigma_b[k] <- sqrt(1/tau_b[k]) # sd for random intercept
    tau[k] ~ dgamma(1,0.1)
    sigma[k] <- sqrt(1/tau[k])
  }
  
    a ~ dnorm(0,0.1)
    base_haz~ dgamma(1,0.1)
  
  for (i in 1:n){b[i] ~ dnorm(0,tau_b[zeta[i]])}
  
  for (i in 1:N){
    mu[i]<-a*t[i]+b[id[i]]+intercept[zeta[id[i]]]
    y[i] ~ dnorm(mu[i],tau[zeta[id[i]]])
  }
  
  for (i in 1:n){
    y_estimated[i]<-Surv_t[i]*a+b[id_Surv[i]]+intercept[zeta[id_Surv[i]]]
    
    elinpred[i] <- exp(beta[zeta[id_Surv[i]]]*x[i]+m[zeta[id_Surv[i]]]*y_estimated[i])
    
    logHaz[i] <-log(base_haz*(Surv_t[i]^(base_haz-1))*elinpred[i])
    
    ### integral of the survival function
    for(l in 1:N_grid){
      f[i,l]<-base_haz*(((Surv_t_grid_point[i,l]+Surv_t_grid_point[i,l+1])/2)^(base_haz-1))*
        exp(beta[zeta[id_Surv[i]]]*x[i]+
              m[zeta[id_Surv[i]]]*
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
",file="DP_weibull_joint.txt")

K=4
N_grid=120

Surv_t_grid_point<-matrix(0,nrow =nrow(dat_surv), ncol = N_grid+1)
for (i in 1:nrow(dat_surv)) {
  Surv_t_grid_point[i,1:(N_grid+1)]<-seq(from=0, to=dat_surv$event_time[i], 
                                         by=dat_surv$event_time[i]/N_grid)
}
n=nrow(dat_surv)
data_test<-list(K=K,
                n=nrow(dat_surv),
                #lmm
                id=dat_lmm$id,
                N=nrow(dat_lmm),
                t=dat_lmm$t,
                y=dat_lmm$y,
                #surv
                Surv_t=dat_surv$event_time,
                id_Surv=dat_surv$id,
                x=dat_surv$x,
                zeros=rep(0,nrow(dat_surv)),
                delta=dat_surv$status,
                #integral
                N_grid=N_grid,
                Surv_t_grid_point=Surv_t_grid_point)

set.seed(SEED.START*1000+SEED.START^2-SEED.START)
jags.m <- jags( model.file= "DP_weibull_joint.txt", 
                data=data_test, 
                inits=list(list(rho=rgamma(1,1,1),
                                beta=rnorm(K,0,1),
                                m=rnorm(K,0,1),base_haz=rgamma(1,1,1),
                                tau_b=rgamma(K,1,0.1),a=rnorm(1,0,1),
                                tau=rgamma(K,1,0.1),intercept=rnorm(K,0,1))),
                n.chains=1, 
                n.iter = 10000,
                n.burnin =8000,
                n.thin=1,
                parameters.to.save=c("base_haz","beta","m","pi",
                                     "sigma_b","sigma","intercept","a","zeta")
)

save(jags.m, file = paste0("DP_all_sim_Comp",SEED.START,".rdata"))


