library(survival)
library(tidyverse)
library(R2jags)
#build function
pred_surv_fun<-function(T_event,sigma_b,a,alpha,beta,m,N=1000,x,intercept){
  b_n<-rnorm(N,0,sigma_b)
  h_t_dependent<- function(s) {alpha * s^(alpha - 1) * exp(a*m*s)}
  H_t_dependent<- integrate(h_t_dependent, lower = 0, upper = T_event)$value
  S_T_b_n<-exp(-exp(beta*x+m*b_n+m*intercept)*H_t_dependent)
  S_T<-mean(S_T_b_n)
  S_T
}
set.seed(111)
#######HDP
xi.1<-rep(0,100)
xi.2<-rep(0,100)

for (j in 1:100) {
load(paste0("Hier_weibull_sim_Comp",j,".rdata"))
load(paste0("dat_train",j,".rdata"))
#data 
dat_surv<-dat%>%group_by(id)%>%slice(1)%>%ungroup()%>%
  select(-y,-t)
dat_surv.1<-dat_surv%>%filter(J==1)
dat_surv.2<-dat_surv%>%filter(J==2)
#KM
fit.1 <- survfit(Surv(event_time, status) ~ 1, data = dat_surv.1)
fit.2 <- survfit(Surv(event_time, status) ~ 1, data = dat_surv.2)
KM.1<-summary(fit.1, times = dat_surv.1$event_time)$surv
KM.2<-summary(fit.2, times = dat_surv.2$event_time)$surv
#survival prob from model
sigma_b<-as.numeric(jags.m$BUGSoutput$mean$sigma_b)
a<-as.numeric(jags.m$BUGSoutput$mean$a)
alpha<-as.numeric(jags.m$BUGSoutput$mean$base_haz)
beta<-as.numeric(jags.m$BUGSoutput$mean$beta)
m<-as.numeric(jags.m$BUGSoutput$mean$m)
intercept<-as.numeric(jags.m$BUGSoutput$mean$intercept)
S.1<-c()
S.2<-c()
pi.1<-jags.m$BUGSoutput$mean$pi[1,]
pi.2<-jags.m$BUGSoutput$mean$pi[2,]
for (i in 1:nrow(dat_surv.1)) {
  S.1[i]<-
    pi.1[1]*pred_surv_fun(T_event =dat_surv.1$event_time[i],
                        sigma_b=sigma_b[1],a=a,
                        alpha=alpha,beta=beta[1],m=m[1],N=1000,
                        x=dat_surv.1$x[i],intercept=intercept[1])+
    pi.1[2]*pred_surv_fun(T_event =dat_surv.1$event_time[i],
                          sigma_b=sigma_b[2],a=a,
                          alpha=alpha,beta=beta[2],m=m[2],N=1000,
                          x=dat_surv.1$x[i],intercept=intercept[2])+
    pi.1[3]*pred_surv_fun(T_event =dat_surv.1$event_time[i],
                          sigma_b=sigma_b[3],a=a,
                          alpha=alpha,beta=beta[3],m=m[3],N=1000,
                          x=dat_surv.1$x[i],intercept=intercept[3])+
    pi.1[4]*pred_surv_fun(T_event =dat_surv.1$event_time[i],
                          sigma_b=sigma_b[4],a=a,
                          alpha=alpha,beta=beta[4],m=m[4],N=1000,
                          x=dat_surv.1$x[i],intercept=intercept[4])
}
for (i in 1:nrow(dat_surv.2)) {
  S.2[i]<-
    pi.2[1]*pred_surv_fun(T_event =dat_surv.2$event_time[i],
                          sigma_b=sigma_b[1],a=a,
                          alpha=alpha,beta=beta[1],m=m[1],N=1000,
                          x=dat_surv.2$x[i],intercept=intercept[1])+
    pi.2[2]*pred_surv_fun(T_event =dat_surv.2$event_time[i],
                          sigma_b=sigma_b[2],a=a,
                          alpha=alpha,beta=beta[2],m=m[2],N=1000,
                          x=dat_surv.2$x[i],intercept=intercept[2])+
    pi.2[3]*pred_surv_fun(T_event =dat_surv.2$event_time[i],
                          sigma_b=sigma_b[3],a=a,
                          alpha=alpha,beta=beta[3],m=m[3],N=1000,
                          x=dat_surv.2$x[i],intercept=intercept[3])+
    pi.2[4]*pred_surv_fun(T_event =dat_surv.2$event_time[i],
                          sigma_b=sigma_b[4],a=a,
                          alpha=alpha,beta=beta[4],m=m[4],N=1000,
                          x=dat_surv.2$x[i],intercept=intercept[4])
}

xi.1[j]<-mean(abs(S.1-KM.1))
xi.2[j]<-mean(abs(S.2-KM.2))
}

gof.hdp<-data.frame(cbind(xi.1,xi.2,
                          xi.all=(xi.1+xi.2)/2))
########DP_ALL
xi.1<-rep(0,100)
xi.2<-rep(0,100)

for (j in 1:100) {
  load(paste0("DP_all_sim_Comp",j,".rdata"))
  load(paste0("/Users/bohanjia/Desktop/research/gsr/Sim_Final/sim_100_train/dat_train",j,".rdata"))
  #data 
  dat_surv<-dat%>%group_by(id)%>%slice(1)%>%ungroup()%>%
    select(-y,-t)
  dat_surv.1<-dat_surv%>%filter(J==1)
  dat_surv.2<-dat_surv%>%filter(J==2)
  #KM
  fit.1 <- survfit(Surv(event_time, status) ~ 1, data = dat_surv.1)
  fit.2 <- survfit(Surv(event_time, status) ~ 1, data = dat_surv.2)
  KM.1<-summary(fit.1, times = dat_surv.1$event_time)$surv
  KM.2<-summary(fit.2, times = dat_surv.2$event_time)$surv
  #survival prob from model
  sigma_b<-as.numeric(jags.m$BUGSoutput$mean$sigma_b)
  a<-as.numeric(jags.m$BUGSoutput$mean$a)
  alpha<-as.numeric(jags.m$BUGSoutput$mean$base_haz)
  beta<-as.numeric(jags.m$BUGSoutput$mean$beta)
  m<-as.numeric(jags.m$BUGSoutput$mean$m)
  intercept<-as.numeric(jags.m$BUGSoutput$mean$intercept)
  S.1<-c()
  S.2<-c()
  pi<-jags.m$BUGSoutput$mean$pi
  for (i in 1:nrow(dat_surv.1)) {
    S.1[i]<-
      pi[1]*pred_surv_fun(T_event =dat_surv.1$event_time[i],
                            sigma_b=sigma_b[1],a=a,
                            alpha=alpha,beta=beta[1],m=m[1],N=1000,
                            x=dat_surv.1$x[i],intercept=intercept[1])+
      pi[2]*pred_surv_fun(T_event =dat_surv.1$event_time[i],
                            sigma_b=sigma_b[2],a=a,
                            alpha=alpha,beta=beta[2],m=m[2],N=1000,
                            x=dat_surv.1$x[i],intercept=intercept[2])+
      pi[3]*pred_surv_fun(T_event =dat_surv.1$event_time[i],
                            sigma_b=sigma_b[3],a=a,
                            alpha=alpha,beta=beta[3],m=m[3],N=1000,
                            x=dat_surv.1$x[i],intercept=intercept[3])+
      pi[4]*pred_surv_fun(T_event =dat_surv.1$event_time[i],
                            sigma_b=sigma_b[4],a=a,
                            alpha=alpha,beta=beta[4],m=m[4],N=1000,
                            x=dat_surv.1$x[i],intercept=intercept[4])
  }
  for (i in 1:nrow(dat_surv.2)) {
    S.2[i]<-
      pi[1]*pred_surv_fun(T_event =dat_surv.2$event_time[i],
                            sigma_b=sigma_b[1],a=a,
                            alpha=alpha,beta=beta[1],m=m[1],N=1000,
                            x=dat_surv.2$x[i],intercept=intercept[1])+
      pi[2]*pred_surv_fun(T_event =dat_surv.2$event_time[i],
                            sigma_b=sigma_b[2],a=a,
                            alpha=alpha,beta=beta[2],m=m[2],N=1000,
                            x=dat_surv.2$x[i],intercept=intercept[2])+
      pi[3]*pred_surv_fun(T_event =dat_surv.2$event_time[i],
                            sigma_b=sigma_b[3],a=a,
                            alpha=alpha,beta=beta[3],m=m[3],N=1000,
                            x=dat_surv.2$x[i],intercept=intercept[3])+
      pi[4]*pred_surv_fun(T_event =dat_surv.2$event_time[i],
                            sigma_b=sigma_b[4],a=a,
                            alpha=alpha,beta=beta[4],m=m[4],N=1000,
                            x=dat_surv.2$x[i],intercept=intercept[4])
  }
  
  xi.1[j]<-mean(abs(S.1-KM.1))
  xi.2[j]<-mean(abs(S.2-KM.2))
}

gof.dp.all<-data.frame(cbind(xi.1,xi.2,
                          xi.all=(xi.1+xi.2)/2))
#######DP_EACH
xi.1<-rep(0,100)
xi.2<-rep(0,100)

for (j in 1:100) {
  load(paste0("/Users/bohanjia/Desktop/research/gsr/Sim_Final/sim_100_train/dat_train",j,".rdata"))
  #data 
  dat_surv<-dat%>%group_by(id)%>%slice(1)%>%ungroup()%>%
    select(-y,-t)
  dat_surv.1<-dat_surv%>%filter(J==1)
  dat_surv.2<-dat_surv%>%filter(J==2)
  #KM
  fit.1 <- survfit(Surv(event_time, status) ~ 1, data = dat_surv.1)
  fit.2 <- survfit(Surv(event_time, status) ~ 1, data = dat_surv.2)
  KM.1<-summary(fit.1, times = dat_surv.1$event_time)$surv
  KM.2<-summary(fit.2, times = dat_surv.2$event_time)$surv
  #survival prob from model
  load(paste0("DP_Each1_sim_Comp",j,".rdata"))
  sigma_b<-as.numeric(jags.m$BUGSoutput$mean$sigma_b)
  a<-as.numeric(jags.m$BUGSoutput$mean$a)
  alpha<-as.numeric(jags.m$BUGSoutput$mean$base_haz)
  beta<-as.numeric(jags.m$BUGSoutput$mean$beta)
  m<-as.numeric(jags.m$BUGSoutput$mean$m)
  intercept<-as.numeric(jags.m$BUGSoutput$mean$intercept)
  pi<-jags.m$BUGSoutput$mean$pi
  S.1<-c()
  for (i in 1:nrow(dat_surv.1)) {
    S.1[i]<-
      pi[1]*pred_surv_fun(T_event =dat_surv.1$event_time[i],
                          sigma_b=sigma_b[1],a=a,
                          alpha=alpha,beta=beta[1],m=m[1],N=1000,
                          x=dat_surv.1$x[i],intercept=intercept[1])+
      pi[2]*pred_surv_fun(T_event =dat_surv.1$event_time[i],
                          sigma_b=sigma_b[2],a=a,
                          alpha=alpha,beta=beta[2],m=m[2],N=1000,
                          x=dat_surv.1$x[i],intercept=intercept[2])+
      pi[3]*pred_surv_fun(T_event =dat_surv.1$event_time[i],
                          sigma_b=sigma_b[3],a=a,
                          alpha=alpha,beta=beta[3],m=m[3],N=1000,
                          x=dat_surv.1$x[i],intercept=intercept[3])+
      pi[4]*pred_surv_fun(T_event =dat_surv.1$event_time[i],
                          sigma_b=sigma_b[4],a=a,
                          alpha=alpha,beta=beta[4],m=m[4],N=1000,
                          x=dat_surv.1$x[i],intercept=intercept[4])
  }
  load(paste0("DP_Each2_sim_Comp",j,".rdata"))
  sigma_b<-as.numeric(jags.m$BUGSoutput$mean$sigma_b)
  a<-as.numeric(jags.m$BUGSoutput$mean$a)
  alpha<-as.numeric(jags.m$BUGSoutput$mean$base_haz)
  beta<-as.numeric(jags.m$BUGSoutput$mean$beta)
  m<-as.numeric(jags.m$BUGSoutput$mean$m)
  intercept<-as.numeric(jags.m$BUGSoutput$mean$intercept)
  pi<-jags.m$BUGSoutput$mean$pi
  S.2<-c()
  for (i in 1:nrow(dat_surv.2)) {
    S.2[i]<-
      pi[1]*pred_surv_fun(T_event =dat_surv.2$event_time[i],
                          sigma_b=sigma_b[1],a=a,
                          alpha=alpha,beta=beta[1],m=m[1],N=1000,
                          x=dat_surv.2$x[i],intercept=intercept[1])+
      pi[2]*pred_surv_fun(T_event =dat_surv.2$event_time[i],
                          sigma_b=sigma_b[2],a=a,
                          alpha=alpha,beta=beta[2],m=m[2],N=1000,
                          x=dat_surv.2$x[i],intercept=intercept[2])+
      pi[3]*pred_surv_fun(T_event =dat_surv.2$event_time[i],
                          sigma_b=sigma_b[3],a=a,
                          alpha=alpha,beta=beta[3],m=m[3],N=1000,
                          x=dat_surv.2$x[i],intercept=intercept[3])+
      pi[4]*pred_surv_fun(T_event =dat_surv.2$event_time[i],
                          sigma_b=sigma_b[4],a=a,
                          alpha=alpha,beta=beta[4],m=m[4],N=1000,
                          x=dat_surv.2$x[i],intercept=intercept[4])
  }
  
  xi.1[j]<-mean(abs(S.1-KM.1))
  xi.2[j]<-mean(abs(S.2-KM.2))
}

gof.dp.each<-data.frame(cbind(xi.1,xi.2,
                         xi.all=(xi.1+xi.2)/2))
save(gof.hdp, file = paste0("gof_hdp.rdata"))
save(gof.dp.each, file = paste0("gof_dp_each.rdata"))
save(gof.dp.all, file = paste0("gof_dp_all.rdata"))

load(paste0("gof_hdp.rdata"))
load(paste0("gof_dp_each.rdata"))
load(paste0("gof_dp_all.rdata"))

colMeans(gof.hdp)
colMeans(gof.dp.each)
colMeans(gof.dp.all)

diff1<-gof.hdp$xi.all-gof.dp.all$xi.all
diff2<-gof.hdp$xi.all-gof.dp.each$xi.all
sd(gof.hdp$xi.all)
sd(gof.dp.each$xi.all)
sd(gof.dp.all$xi.all)
#whole
wilcox.test(gof.hdp$xi.all,gof.dp.each$xi.all,
            paired = TRUE, exact = FALSE) 

wilcox.test(gof.hdp$xi.all, gof.dp.all$xi.all,
            paired = TRUE, exact = FALSE)
