SEED.START = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
library(R2jags)
library(tidyverse)
load(paste0("dat_CD4_combine.rdata"))
load(paste0("dat_surv_combine.rdata"))
dat_CD4_combine<-dat_CD4_combine%>%
  select(-pop)
dat_CD4_combine<-inner_join(dat_CD4_combine,
                            dat_surv_combine,by="id")
dat_CD4_combine<-dat_CD4_combine%>%
  mutate(post_event=visit<ART_init)
dat_CD4_combine$post_event<-as.numeric(dat_CD4_combine$post_event)
# model
cat("model{
  
  for (k in 1:K)
  {
    ###survival
    m[k] ~ dnorm(0,0.1) # parameters for y
    #base_haz[k] ~ dgamma(0.1,0.1)
    
    ###lmm
    intercept[k] ~ dnorm(0,0.1) # fixed intercept
    tau_b[k] ~ dgamma(0.1,0.1)
    sigma_b[k] <- sqrt(1/tau_b[k]) # sd for random intercept
  }
  
  base_haz ~ dgamma(0.1,0.1)
  tau ~ dgamma(0.1,0.1)
  sigma <- sqrt(1/tau)
  a[1] ~ dnorm(0,0.1) # fixed effect for t
  a[2] ~ dnorm(0,0.1) # fixed effect for t
  
  for (i in 1:n){b[i] ~ dnorm(0,tau_b[zeta[i]])}
  
  for (i in 1:N){
    
    mu[i]<-(a[1]+a[2]*post_event[i])*t[i]+b[id[i]]+intercept[zeta[id[i]]]
    y[i] ~ dnorm(mu[i],tau)
  }
  
  for (i in 1:N_Surv){
    y_estimated[i]<-Surv_t[i]*a[1]+b[id_Surv[i]]+intercept[zeta[id_Surv[i]]]
    
    elinpred[i] <- exp(m[zeta[id_Surv[i]]]*y_estimated[i])
    
    logHaz[i] <-log(base_haz*(Surv_t[i]^(base_haz-1))*elinpred[i])
    
    ### integral of the survival function
    for(l in 1:N_grid){
      f[i,l]<-base_haz*(((Surv_t_grid_point[i,l]+Surv_t_grid_point[i,l+1])/2)^(base_haz-1))*
        exp(m[zeta[id_Surv[i]]]*
              (b[id_Surv[i]]+intercept[zeta[id_Surv[i]]]))*
        exp(m[zeta[id_Surv[i]]]*a[1]*((Surv_t_grid_point[i,l]+Surv_t_grid_point[i,l+1])/2))*
        (Surv_t_grid_point[i,l+1]-Surv_t_grid_point[i,l])}
    
    logSurv[i] <- -sum(f[i,1:N_grid])
    
    phi[i]<- 100 - delta[i]*logHaz[i]-logSurv[i]
    zeros[i]~dpois(phi[i])
  }
  
  for (i in 1:n_k){zeta[i] ~ dcat(pi[1,])}
  for (i in (n_k+1):n){zeta[i] ~ dcat(pi[2,])}
  
  rho ~ dgamma(2,2)
  for (k in 1:(K-1)) { V[k] ~ dbeta(1,rho) T(0.001,0.999) }
  B[1] <- V[1]
  R[1] <- 1 - V[1]
  
  for (k in 2:(K-1)) {
  B[k] <- V[k] * R[k-1]
  R[k] <- R[k-1] * (1 - V[k])}
  
  B[K] <- R[K-1]
  
  alpha ~ dgamma(2,2)
  for (j in 1:J){
  for (k in 1:(K-1)){
    shape1[j,k] <- max(alpha*B[k],0.001)
    shape2[j,k] <- max(alpha*sum(B[(k+1):K]), 0.001)
    U[j,k] ~ dbeta(shape1[j,k], shape2[j,k]) T(0.001,0.999)
  }
  U[j,K] <- 1
}
  
  for (j in 1:J){
    pi[j,1]<-U[j,1]
    for (k in 2:K){pi[j,k]<-U[j,k]*(1-U[j,k-1])*pi[j,k-1]/U[j,k-1]}
  }
  
  
  
}
",file="HDP_weibull_post_joint.txt")

K=10
J=2
N_grid=100

Surv_t_grid_point<-matrix(0,nrow =nrow(dat_surv_combine), ncol = N_grid+1)
for (i in 1:nrow(dat_surv_combine)) {
  Surv_t_grid_point[i,1:(N_grid+1)]<-seq(from=0, to=dat_surv_combine$ART_init[i], 
                                         by=dat_surv_combine$ART_init[i]/N_grid)
}
n=nrow(dat_surv_combine)
n_594<-sum(dat_surv_combine$pop=="594")
data_test<-list(K=K,
                n=nrow(dat_surv_combine),
                n_k=n_594,
                J=J,
                #lmm
                id=dat_CD4_combine$id,
                N=nrow(dat_CD4_combine),
                t=dat_CD4_combine$visit,
                y=dat_CD4_combine$cd4count,
                post_event=dat_CD4_combine$post_event,
                #surv
                Surv_t=dat_surv_combine$ART_init,
                id_Surv=dat_surv_combine$id,
                N_Surv=nrow(dat_surv_combine),
                zeros=rep(0,nrow(dat_surv_combine)),
                delta=dat_surv_combine$status,
                #integral
                N_grid=N_grid,
                Surv_t_grid_point=Surv_t_grid_point)

set.seed(SEED.START^2+200*SEED.START)
jags.m <- jags( model.file= "HDP_weibull_post_joint.txt", 
                data=data_test, 
                inits=list(list(rho=rgamma(1,1,1),alpha=rgamma(1,1,1),
                                m=rnorm(K,0,1),base_haz=rgamma(1,1,0.5),
                                tau_b=rgamma(K,1,0.1),a=rnorm(2,0,1),
                                tau=rgamma(1,1,0.1),intercept=rnorm(K,0,1))),
                n.chains=1, 
                n.iter = 20000,
                n.burnin =17000,
                n.thin=1,
                parameters.to.save=c("base_haz","m","pi",
                                     "sigma_b","sigma","intercept","a","zeta")
)
save(jags.m, file = paste0("594320_HDP_Post.rdata"))
