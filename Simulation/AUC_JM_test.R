library(tidyverse)
library(R2jags)
#function
P_u_given_tYk<-function(t,u,Y,t_Y,
                        a,intercept,sigma,sigma_b,
                        x,alpha,beta,m,
                        N=1000){
  b_n<-rnorm(N,0,sigma_b)
  log_P_Y<-c()
  log_S_t<-c()
  log_S_u<-c()
  h_t_dependent<- function(s) {alpha * s^(alpha - 1) * exp(a*m*s)}
  H_t_dependent<- integrate(h_t_dependent, lower = 0, upper = t)$value
  H_u_dependent<- integrate(h_t_dependent, lower = 0, upper = u)$value
  #if we have random intercept then we need to put it in the loop
  for (n in 1:N) {
    log_P_Y[n]<-sum(dnorm(Y,t_Y*a+b_n[n]+intercept,sigma,log = TRUE))
    log_S_t[n]<- -exp(beta*x+m*b_n[n]+m*intercept)*H_t_dependent
    log_S_u[n]<- -exp(beta*x+m*b_n[n]+m*intercept)*H_u_dependent}
  
  logsumexp <- function(z) {
    zmax <- max(z)
    zmax + log(sum(exp(z - zmax)))
  }
  
  log_num_terms <- log_P_Y + log_S_u
  log_den_terms <- log_P_Y + log_S_t
  
  log_num <- logsumexp(log_num_terms)
  log_den <- logsumexp(log_den_terms)
  
  exp(log_num - log_den)
}
# t and u
t_given<-c(2,6,12,18)
delta_t<-c(1,2,3)
AUC_JM_sim<-list()
for (s in 1:100) {
  SN_SP<-data.frame(C=seq(from=0,to=1,by=0.001))
  load(paste0("JM_sim_Comp",s,".rdata"))
  load(paste0("dat_test",s,".rdata"))
  #para
  sigma_b<-as.numeric(jags.m$BUGSoutput$mean$sigma_b)
  sigma<-as.numeric(jags.m$BUGSoutput$mean$sigma)
  a<-as.numeric(jags.m$BUGSoutput$mean$a)
  alpha<-as.numeric(jags.m$BUGSoutput$mean$base_haz)
  beta<-as.numeric(jags.m$BUGSoutput$mean$beta)
  m<-as.numeric(jags.m$BUGSoutput$mean$m)
  intercept<-as.numeric(jags.m$BUGSoutput$mean$intercept)
  #calculation
  for (j in 1:length(t_given)) {
    #data
    dat_lmm<-dat%>%filter(event_time>=t)%>%
      filter(t<=t_given[j])%>%select(-event_time,-status)
    dat_surv<-dat%>%filter(event_time>=t)%>%
      filter(t<=t_given[j])%>%group_by(id)%>%
      slice(1)%>%ungroup()%>%
      select(-y,-t)
    id.left<-dat_surv$id
    for (l in 1:length(delta_t)) {
      marg_P_u_given_tY<-rep(0,nrow(dat_surv))
      for (i in 1:nrow(dat_surv)) {
        Y_i<-data.frame(dat_lmm%>%filter(id==id.left[i]))
        marg_P_u_given_tY[i]<-P_u_given_tYk(t=t_given[j],u=t_given[j]+delta_t[l],Y=Y_i$y,t_Y = Y_i$t,
                                            a=a,intercept=intercept,sigma=sigma,sigma_b = sigma_b,
                                            x=dat_surv$x[i],alpha=alpha,beta=beta,m=m,
                                            N=1000)
      } 
      #AUC
      C<-seq(from=0,to=1,by=0.001)
      I_survivor<-dat_surv$event_time>t_given[j]
      # "dat_surv$event_time>t_given[j]+delta_t[l]" and 
      # "dat_surv$event_time>(t_given[j]+delta_t[l])" are the same in R
      Omega<-dat_surv$event_time<=t_given[j]+delta_t[l] # we do not have T_I < u and \delta_i=0
      SN<-c()
      for (c in 1:length(C)) {
        SN[c]<- sum(I_survivor*(marg_P_u_given_tY<=C[c])*Omega)/sum(I_survivor*Omega)
      }
      Phi<-dat_surv$event_time>t_given[j]+delta_t[l] # we do not have T_I < u and \delta_i=0
      SP<-c()
      for (c in 1:length(C)) {
        SP[c]<- sum(I_survivor*(marg_P_u_given_tY>C[c])*Phi)/sum(I_survivor*Phi)
      }
      SN_SP<-cbind(SN_SP,data.frame(SN,SP))
      colnames(SN_SP)[c(ncol(SN_SP)-1,ncol(SN_SP))]<-
        c(paste0("SN-",t_given[j],t_given[j]+delta_t[l]), 
          paste0("SP-",t_given[j],t_given[j]+delta_t[l]))
    }}
  rownames(SN_SP)<-1:nrow(SN_SP)
  AUC_JM_sim[[s]]<-SN_SP
}
save(AUC_JM_sim,file = paste0("AUC_JM_sim.rdata"))

