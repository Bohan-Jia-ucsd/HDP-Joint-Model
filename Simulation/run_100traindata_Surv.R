
SEED.START = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(SEED.START*100+SEED.START^2+SEED.START)
library(tidyverse)
library(R2jags)
load( paste0("dat_train",SEED.START,".rdata"))
dat_lmm<-dat%>%select(-event_time,-status)
dat_surv<-dat%>%group_by(id)%>%slice(1)%>%ungroup()%>%
  select(-y,-t)

dat_combine<-merge(dat_surv[,c(-2,-5,-6)],dat_lmm,by="id")
dat_combine<-dat_combine%>%filter(t<=event_time)
n_item<-dat_combine%>%group_by(id)%>%summarise(n=n())
id_new<-c()
for (i in 1:nrow(n_item)) {
  id_new<-c(id_new,rep(i,n_item$n[i]))
}
dat_combine$id<-id_new
n_aid<-length(unique(dat_combine[dat_combine$J==1,"id"]))
dat_surv<-dat_combine%>%group_by(id)%>%slice(n())%>%ungroup()
s_mat<-matrix(0,nrow = nrow(dat_surv),ncol = 5+1)
y_mat<-matrix(0,nrow = nrow(dat_surv),ncol = 5+1)
for (i in 1:nrow(dat_surv)) {
  for (j in 1:n_item$n[i]) {
    s_mat[i,j]<-dat_combine[dat_combine$id==i,"t"][j]
    y_mat[i,j]<-dat_combine[dat_combine$id==i,"y"][j]
  }
  s_mat[i,n_item$n[i]+1]<-dat_surv$event_time[i]
}

# model
cat("model{
  
  for (k in 1:K)
  {
    ###survival
    m[k] ~ dnorm(0,0.1) # parameters for y
    beta[k] ~ dnorm(0,0.1) # parameters for baseline covariates
  }
  
  base_haz ~ dgamma(0.1,0.1)
  
  for (i in 1:N_Surv){
    
    elinpred[i] <- exp(beta[zeta[id_Surv[i]]]*x[i]+m[zeta[id_Surv[i]]]*x_T[i])
    
    logHaz[i] <-log(base_haz*(Surv_t[i]^(base_haz-1))*elinpred[i])
    
    ### integral of the survival function
    for(l in 1:N_int[i]){
      f[i,l]<-
      (s[i,l+1]^(base_haz)-s[i,l]^(base_haz))*
      exp(beta[zeta[id_Surv[i]]]*x[i]+m[zeta[id_Surv[i]]]*y[i,l])
      }
    
    logSurv[i] <- -sum(f[i,1:N_int[i]])
    
    phi[i]<- 100- delta[i]*logHaz[i]-logSurv[i]
    zeros[i]~dpois(phi[i])
  }
  
  for (i in 1:n_k){zeta[i] ~ dcat(pi[1,])}
  for (i in (n_k+1):n){zeta[i] ~ dcat(pi[2,])}
  
  rho ~ dgamma(1,1)
  for (k in 1:(K-1)) { V[k] ~ dbeta(1,rho) T(0.0001,0.9999) }
  V[K] <- 1
  B[1] <- V[1]
  
  for (k in 2:K) {
    B[k] <- V[k] * (1-V[k-1]) * B[k-1]/V[k-1]
  }
  
  alpha ~ dgamma(1,1)
  for (j in 1:J){
    for (k in 1:(K-1)){
      U[j,k] ~ dbeta(alpha*B[k],alpha*sum(B[(k+1):K])) T(0.001,0.999)
    }
    U[j,K]<-1
  }
  
  for (j in 1:J){
    pi[j,1]<-U[j,1]
    for (k in 2:K){pi[j,k]<-U[j,k]*(1-U[j,k-1])*pi[j,k-1]/U[j,k-1]}
  }
  
  
  
}
",file="HDP_Surv_joint.txt")

K=4
J=2

data_test<-list(K=K,
                n=nrow(dat_surv),
                n_k=n_aid,
                J=J,
                #surv
                x_T=dat_surv$y,
                N_int=n_item$n,
                s=s_mat,
                y=y_mat,
                Surv_t=dat_surv$event_time,
                id_Surv=dat_surv$id,
                N_Surv=nrow(dat_surv),
                x=dat_surv$x,
                zeros=rep(0,nrow(dat_surv)),
                delta=dat_surv$status)


jags.m <- jags( model.file= "HDP_Surv_joint.txt", 
                data=data_test, 
                inits=list(list(rho=rgamma(1,1,1),alpha=rgamma(1,1,1),
                                beta=rnorm(K,0,1),
                                m=rnorm(K,0,1),base_haz=rgamma(1,1,0.5))),
                n.chains=1, 
                n.iter = 10000,
                n.burnin =8000,
                n.thin=1,
                parameters.to.save=c("base_haz","beta","m","pi","zeta")
)
save(jags.m, file = paste0("Surv_HDP_sim_Comp",SEED.START,".rdata"))
