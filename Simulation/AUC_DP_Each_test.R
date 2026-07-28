library(tidyverse)
library(R2jags)
#function
# 1
pi_k_given_tY_propto<-function(t,Y,t_Y,pi,
                               a,intercept,sigma,sigma_b,
                               x,alpha,beta,m,
                               N=1000){
  b_n<-rnorm(N,0,sigma_b)
  P_Y<-c()
  S_t<-c()
  for (n in 1:N) {
    P_Y[n]<-prod(dnorm(Y,t_Y*a+b_n[n]+intercept,sigma))
    h_t_dependent<- function(s) {alpha * s^(alpha - 1) * exp(a*m*s)}
    H_t_dependent<- integrate(h_t_dependent, lower = 0, upper = t)$value
    S_t[n]<-exp(-exp(beta*x+m*b_n[n]+m*intercept)*H_t_dependent)
  }
  mean(P_Y*S_t)*pi
}
# 2
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
AUC_DP_Each_sim<-list()
for (s in 1:100) {
  SN_SP<-data.frame(C=seq(from=0,to=1,by=0.001))
  load(paste0("dat_test",s,".rdata"))
  #para
  load(paste0("DP_Each1_sim_Comp",s,".rdata"))
  sigma_b<-as.numeric(jags.m$BUGSoutput$mean$sigma_b)
  sigma<-as.numeric(jags.m$BUGSoutput$mean$sigma)
  a<-as.numeric(jags.m$BUGSoutput$mean$a)
  alpha<-as.numeric(jags.m$BUGSoutput$mean$base_haz)
  beta<-as.numeric(jags.m$BUGSoutput$mean$beta)
  m<-as.numeric(jags.m$BUGSoutput$mean$m)
  intercept<-as.numeric(jags.m$BUGSoutput$mean$intercept)
  pi<-as.numeric(jags.m$BUGSoutput$mean$pi)
  load(paste0("DP_Each2_sim_Comp",s,".rdata"))
  sigma_b<-rbind(sigma_b,as.numeric(jags.m$BUGSoutput$mean$sigma_b))
  sigma<-rbind(sigma,as.numeric(jags.m$BUGSoutput$mean$sigma))
  a<-c(a,as.numeric(jags.m$BUGSoutput$mean$a))
  alpha<-c(alpha,as.numeric(jags.m$BUGSoutput$mean$base_haz))
  beta<-rbind(beta,as.numeric(jags.m$BUGSoutput$mean$beta))
  m<-rbind(m,as.numeric(jags.m$BUGSoutput$mean$m))
  intercept<-rbind(intercept,as.numeric(jags.m$BUGSoutput$mean$intercept))
  pi<-rbind(pi,as.numeric(jags.m$BUGSoutput$mean$pi))
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
        pi_k_given_tY<-rep(0,4)
        pi_k_given_tY[1]<-pi_k_given_tY_propto(t=t_given[j],Y=Y_i$y,t_Y = Y_i$t,
                                               pi=pi[dat_surv$J[i],1],a=a[dat_surv$J[i]],
                                               intercept=intercept[dat_surv$J[i],1],
                                               sigma=sigma[dat_surv$J[i],1],sigma_b = sigma_b[dat_surv$J[i],1],
                                               x=dat_surv$x[i],alpha=alpha[dat_surv$J[i]],
                                               beta=beta[dat_surv$J[i],1],
                                               m=m[dat_surv$J[i],1],N=1000)
        pi_k_given_tY[2]<-pi_k_given_tY_propto(t=t_given[j],Y=Y_i$y,t_Y = Y_i$t,
                                               pi=pi[dat_surv$J[i],2],a=a[dat_surv$J[i]],
                                               intercept=intercept[dat_surv$J[i],2],
                                               sigma=sigma[dat_surv$J[i],2],sigma_b = sigma_b[dat_surv$J[i],2],
                                               x=dat_surv$x[i],alpha=alpha[dat_surv$J[i]],
                                               beta=beta[dat_surv$J[i],2],
                                               m=m[2],N=1000)
        pi_k_given_tY[3]<-pi_k_given_tY_propto(t=t_given[j],Y=Y_i$y,t_Y = Y_i$t,
                                               pi=pi[dat_surv$J[i],3],a=a[dat_surv$J[i]],
                                               intercept=intercept[dat_surv$J[i],3],
                                               sigma=sigma[dat_surv$J[i],3],sigma_b = sigma_b[dat_surv$J[i],3],
                                               x=dat_surv$x[i],alpha=alpha[dat_surv$J[i]],
                                               beta=beta[dat_surv$J[i],3],
                                               m=m[3],N=1000)
        pi_k_given_tY[4]<-pi_k_given_tY_propto(t=t_given[j],Y=Y_i$y,t_Y = Y_i$t,
                                               pi=pi[dat_surv$J[i],4],a=a[dat_surv$J[i]],
                                               intercept=intercept[dat_surv$J[i],4],
                                               sigma=sigma[dat_surv$J[i],4],sigma_b = sigma_b[dat_surv$J[i],4],
                                               x=dat_surv$x[i],alpha=alpha[dat_surv$J[i]],
                                               beta=beta[dat_surv$J[i],4],
                                               m=m[4],N=1000)
        pi_k_given_tY<-pi_k_given_tY/sum(pi_k_given_tY)
        ##
        P_u_given_tYk_store<-rep(0,4)
        P_u_given_tYk_store[1]<-P_u_given_tYk(t=t_given[j],u=t_given[j]+delta_t[l],Y=Y_i$y,t_Y = Y_i$t,
                                              a=a[dat_surv$J[i]],intercept=intercept[dat_surv$J[i],1],
                                              sigma=sigma[dat_surv$J[i],1],sigma_b = sigma_b[dat_surv$J[i],1],
                                              x=dat_surv$x[i],alpha=alpha[dat_surv$J[i]],
                                              beta=beta[dat_surv$J[i],1],m=m[dat_surv$J[i],1],
                                              N=1000)
        P_u_given_tYk_store[2]<-P_u_given_tYk(t=t_given[j],u=t_given[j]+delta_t[l],Y=Y_i$y,t_Y = Y_i$t,
                                              a=a[dat_surv$J[i]],intercept=intercept[dat_surv$J[i],2],
                                              sigma=sigma[dat_surv$J[i],2],sigma_b = sigma_b[dat_surv$J[i],2],
                                              x=dat_surv$x[i],alpha=alpha[dat_surv$J[i]],
                                              beta=beta[dat_surv$J[i],2],m=m[dat_surv$J[i],2],
                                              N=1000)
        P_u_given_tYk_store[3]<-P_u_given_tYk(t=t_given[j],u=t_given[j]+delta_t[l],Y=Y_i$y,t_Y = Y_i$t,
                                              a=a[dat_surv$J[i]],intercept=intercept[dat_surv$J[i],3],
                                              sigma=sigma[dat_surv$J[i],3],sigma_b = sigma_b[dat_surv$J[i],3],
                                              x=dat_surv$x[i],alpha=alpha[dat_surv$J[i]],
                                              beta=beta[dat_surv$J[i],3],m=m[dat_surv$J[i],3],
                                              N=1000)
        P_u_given_tYk_store[4]<-P_u_given_tYk(t=t_given[j],u=t_given[j]+delta_t[l],Y=Y_i$y,t_Y = Y_i$t,
                                              a=a[dat_surv$J[i]],intercept=intercept[dat_surv$J[i],4],
                                              sigma=sigma[dat_surv$J[i],4],sigma_b = sigma_b[dat_surv$J[i],4],
                                              x=dat_surv$x[i],alpha=alpha[dat_surv$J[i]],
                                              beta=beta[dat_surv$J[i],4],m=m[dat_surv$J[i],4],
                                              N=1000)
        marg_P_u_given_tY[i]<-sum(P_u_given_tYk_store*pi_k_given_tY)
      } 
      #AUC
      C<-seq(from=0,to=1,by=0.001)
      I_survivor<-dat_surv$event_time>t_given[j]
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
  AUC_DP_Each_sim[[s]]<-SN_SP
}
save(AUC_DP_Each_sim,file = paste0("AUC_DP_Each_sim.rdata"))
