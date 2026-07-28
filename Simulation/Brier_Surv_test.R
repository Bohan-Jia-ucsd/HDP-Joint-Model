library(tidyverse)
library(R2jags)
# t and u
t_given<-c(2,6,12,18)
delta_t<-c(1,2,3)
Brier_Surv_sim<-list()
for (s in 1:100) {
Brier_Surv<-matrix(0,ncol = length(t_given),nrow = length(delta_t))
colnames(Brier_Surv)<-c(2,6,12,18)
rownames(Brier_Surv)<-c(1,2,3)
load(paste0("Surv_HDP_sim_Comp",s,".rdata"))
load(paste0("dat_test",s,".rdata"))
#para
alpha<-as.numeric(jags.m$BUGSoutput$mean$base_haz)
beta<-as.numeric(jags.m$BUGSoutput$mean$beta)
m<-as.numeric(jags.m$BUGSoutput$mean$m)
pi<-rbind(as.numeric(jags.m$BUGSoutput$mean$pi[1,]),
          as.numeric(jags.m$BUGSoutput$mean$pi[2,]))
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
  n_item<-dat_lmm%>%group_by(id)%>%summarise(n=n())
  for (l in 1:length(delta_t)) {
marg_P_u_given_tY<-rep(0,nrow(dat_surv))
for (i in 1:nrow(dat_surv)) {
  s_vec<-rep(0,ncol = 5+1) # 0 2 6 12 18 +event_time, at most
  y_vec<-rep(0,ncol = 5+1)
  for (e in 1:n_item$n[i]) {
    s_vec[e]<-dat_lmm[dat_lmm$id==id.left[i],"t"][e]
    y_vec[e]<-dat_lmm[dat_lmm$id==id.left[i],"y"][e]
  }
  s_vec[n_item$n[i]+1]<-t_given[j]
  s_vec_u<-rep(0,ncol = 5+1)
  Surv_t<-rep(0,4) # 4 classes
  Surv_u<-rep(0,4)
  log_surv<-rep(0,4)
  for (e in 1:n_item$n[i]) {
    log_surv<-log_surv+(s_vec[e+1]^alpha-s_vec[e]^alpha)*
      c(exp(dat_surv$x[i]*beta[1]+m[1]*y_vec[e]),
        exp(dat_surv$x[i]*beta[2]+m[2]*y_vec[e]),
        exp(dat_surv$x[i]*beta[3]+m[3]*y_vec[e]),
        exp(dat_surv$x[i]*beta[4]+m[4]*y_vec[e]))
    s_vec_u[e]<-dat_lmm[dat_lmm$id==id.left[i],"t"][e]
  }
  Surv_t<-exp(-log_surv)
  s_vec_u[n_item$n[i]+1]<-t_given[j]+delta_t[l]
  log_surv<-rep(0,4)
  for (e in 1:n_item$n[i]) {
    log_surv<-log_surv+(s_vec_u[e+1]^alpha-s_vec_u[e]^alpha)*
      c(exp(dat_surv$x[i]*beta[1]+m[1]*y_vec[e]),
        exp(dat_surv$x[i]*beta[2]+m[2]*y_vec[e]),
        exp(dat_surv$x[i]*beta[3]+m[3]*y_vec[e]),
        exp(dat_surv$x[i]*beta[4]+m[4]*y_vec[e]))
  }
  Surv_u<-exp(-log_surv)
  numerator<-Surv_u*pi[dat_surv$J[i],]
  denominator<-sum(Surv_t*pi[dat_surv$J[i],])
  marg_P_u_given_tY[i]<-sum(numerator)/denominator
}

R_t<-nrow(dat_surv[dat_surv$event_time>t_given[j],])
I_survivor<-dat_surv$event_time>t_given[j]
I_T_u<-dat_surv$event_time>t_given[j]+delta_t[l]
Brier_Surv[l,j]<-sum(I_survivor*
                      (I_T_u*(1-marg_P_u_given_tY)^2+
                                  (1-I_T_u)*dat_surv$status*(0-marg_P_u_given_tY)^2))/R_t
}}
Brier_Surv_sim[[s]]<-Brier_Surv
}
save(Brier_Surv_sim,file = paste0("Brier_Surv_sim.rdata"))
