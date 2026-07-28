
SEED.START = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(SEED.START*100+SEED.START^2+SEED.START)
library(tidyverse)
library(R2jags)
load( paste0("dat_train",SEED.START,".rdata"))
dat_lmm<-dat%>%select(-event_time,-status)
dat_surv<-dat%>%group_by(id)%>%slice(1)%>%ungroup()%>%
  select(-y,-t)

# model
cat("model{
  

    ###survival
    m ~ dnorm(0,0.1) # parameters for y
    beta ~ dnorm(0,0.1) # parameters for baseline covariates
    
    ###lmm
    intercept ~ dnorm(0,0.1) # fixed intercept
    tau_b ~ dgamma(1,0.1)
    sigma_b <- sqrt(1/tau_b) # sd for random intercept
    tau ~ dgamma(1,0.1)
    sigma <- sqrt(1/tau)
  
  
    a ~ dnorm(0,0.1)
    base_haz~ dgamma(1,0.1)
  
  for (i in 1:n){b[i] ~ dnorm(0,tau_b)}
  
  for (i in 1:N){
    mu[i]<-a*t[i]+b[id[i]]+intercept
    y[i] ~ dnorm(mu[i],tau)
  }
  
  for (i in 1:n){
    y_estimated[i]<-Surv_t[i]*a+b[id_Surv[i]]+intercept
    
    elinpred[i] <- exp(beta*x[i]+m*y_estimated[i])
    
    logHaz[i] <-log(base_haz*(Surv_t[i]^(base_haz-1))*elinpred[i])
    
    ### integral of the survival function
    for(l in 1:N_grid){
      f[i,l]<-base_haz*(((Surv_t_grid_point[i,l]+Surv_t_grid_point[i,l+1])/2)^(base_haz-1))*
        exp(beta*x[i]+
              m*(b[id_Surv[i]]+intercept))*
        exp(m*a*((Surv_t_grid_point[i,l]+Surv_t_grid_point[i,l+1])/2))*
        (Surv_t_grid_point[i,l+1]-Surv_t_grid_point[i,l])}
    
    logSurv[i] <- -sum(f[i,1:N_grid])
    
    phi[i]<-  100 - delta[i]*logHaz[i]-logSurv[i]
    zeros[i]~dpois(phi[i])
  }
  
  
  
}
",file="JM.txt")

N_grid=120

Surv_t_grid_point<-matrix(0,nrow =nrow(dat_surv), ncol = N_grid+1)
for (i in 1:nrow(dat_surv)) {
  Surv_t_grid_point[i,1:(N_grid+1)]<-seq(from=0, to=dat_surv$event_time[i], 
                                         by=dat_surv$event_time[i]/N_grid)
}
n=nrow(dat_surv)
data_test<-list(n=nrow(dat_surv),
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


jags.m <- jags( model.file= "JM.txt", 
                data=data_test, 
                inits=list(list(beta=rnorm(1,0,1),
                                m=rnorm(1,0,1),base_haz=rgamma(1,1,1),
                                tau_b=rgamma(1,1,0.1),a=rnorm(1,0,1),
                                tau=rgamma(1,1,0.1),intercept=rnorm(1,0,1))),
                n.chains=1, 
                n.iter = 10000,
                n.burnin =8000,
                n.thin=1,
                parameters.to.save=c("base_haz","beta","m",
                                     "sigma_b","sigma","intercept","a")
)

save(jags.m, file = paste0("JM_sim_Comp",SEED.START,".rdata"))


