library(tidyverse)
####All parameters are based on the real data and DP model!!!!!#########

######parameters##############
#LMM
Itcpt1<- 4.7 
Itcpt2<- 11.0
#same fixed effect on t
beta_t1.1<- -0.1
beta_t2.1<- -0.1
#same fixed effect after the non-terminal event on t
#we can assume that the pre-event parameter of t is different from
#the post-event parameter of t
#(but we do not do that here!!)
beta_t1.2<- -0.1
beta_t2.2<- -0.1
rand_int_sigma1<-1.5
rand_int_sigma2<-3.9
sigma_1<-1.4
sigma_2<-2.7
t_no_missing<-c(0,2,6,12,18)

#Survival
beta_d1<- 0.6 # parameter for drug
beta_d2<- 0 # no treatment effect
assoct1<- -1.6
assoct2<- -0.6
#same baseline hazard weibull parameter
alpha1<-1.3
alpha2<-1.3

######Newton Method Function##############
newton_bisec <- function(alpha,X,b,beta,C,a,intercept,eta=0,tol = 1e-10,
                         max_it=10000) {
  h_t_dependent <- function(s) {
    alpha * s^(alpha - 1) * exp(a_tilde * s)} # first derivative
  
  H_t_dependent <- function(t){
    integrate(h_t_dependent , lower = 0, upper = t)$value} # function we want to solve
  
  beta_tilde<-beta+C*eta
  a_tilde<-C*a
  exp_constant<-exp(X*beta_tilde+C*b+C*intercept)
  u <- runif(1,0,1)
  # when a_tilde < 0, exp(a_tilde * s) is decreasing
  # v might be >= H (inf) !!!!! 
  # in this case, no solution for v=H(t), loop of 
  # " while (H(thi) < v){thi <- 2 * thi}" will never stop!!!!
  # manually discard such v!
  if (a_tilde<0){
    H_inf <- integrate(h_t_dependent, 0, Inf)$value
    v <-  -log(u)/exp_constant
    while(v >= H_inf){
      u <- runif(1,0,1)
      v <-  -log(u)/exp_constant
    }
  }else{
    v <-  -log(u)/exp_constant
  }
  
  # bracket
  thi <- 1
  while (H_t_dependent(thi) < v){thi <- 2 * thi}
  tlo <- 0
  t  <- min(thi, v^(1/alpha))  # initial guess
  val <- H_t_dependent(t) - v
  it<-0
  while (abs(val) > tol && it < max_it) {
    # try Newton
    t_new <- t - val / h_t_dependent(t)
    # fallback to bisection if out of bracket or not improving
    if (t_new <= tlo || t_new >= thi) {t_new <- 0.5 * (tlo + thi)}else{
      t_new <-t_new}
    val_new <- H_t_dependent(t_new) - v
    # update bracket
    if (val_new > 0) { thi <- t_new } else { tlo <- t_new }
    # accept step if improving, else bisect
    if (abs(val_new) <= abs(val)) { 
      t <- t_new
      val <- val_new }else { 
        t <- 0.5 * (tlo + thi)
        val <- H_t_dependent(t) - v }
    it<-it+1
  }
  if (it == max_it) {NA} else{t}
}

#####lets simulate!!!
for (sim.order in 1:100) {
  seed.test<-sim.order*30+sim.order
  set.seed(seed.test)

n_aid<-round(467/3,0) ### 307 from real data
n_no_aid<-round(467/3,0) ### 160 from real data
######aids group##############
random_intercept1<-rnorm(n_aid,0,rand_int_sigma1)
random_intercept2<-rnorm(n_aid,0,rand_int_sigma2)
## pre-event
y.pre<-data.frame()
for (i in 1:n_aid) {
  y_1.pre<-t_no_missing*beta_t1.1+random_intercept1[i]+ #.1
    rnorm(length(t_no_missing),0,sigma_1)+Itcpt1
  y_2.pre<-t_no_missing*beta_t2.1+random_intercept2[i]+
    rnorm(length(t_no_missing),0,sigma_2)+Itcpt2
  y.pre<-rbind(y.pre,
                cbind(id=rep(i,length(t_no_missing)),y_1.pre,y_2.pre,
               t=t_no_missing))
}
## post-event
y.post<-data.frame()
for (i in 1:n_aid) {
  y_1.post<-t_no_missing*beta_t1.2+random_intercept1[i]+ #.2
    rnorm(length(t_no_missing),0,sigma_1)+Itcpt1
  y_2.post<-t_no_missing*beta_t2.2+random_intercept2[i]+
    rnorm(length(t_no_missing),0,sigma_2)+Itcpt2
  y.post<-rbind(y.post,
                cbind(id=rep(i,length(t_no_missing)),y_1.post,y_2.post,
                t=t_no_missing))
}
## event time
x<-sample(c(0,1),n_aid,replace = TRUE) #baseline covariate: drug
t_1<-c()
for (i in 1:n_aid) {
  t_1[i]<-newton_bisec(alpha=alpha1,X=x[i],b=random_intercept1[i],
                       beta=beta_d1,C=assoct1,a=beta_t1.1,intercept=Itcpt1,
                       eta=0,tol = 1e-10)      #.1
}
quantile(t_1,na.rm = T)
t_2<-c()
for (i in 1:n_aid) {
  t_2[i]<-newton_bisec(alpha=alpha2,X=x[i],b=random_intercept2[i],
                       beta=beta_d2,C=assoct2,a=beta_t2.1,intercept=Itcpt2,
                       eta=0,tol = 1e-10)      #.1
}
quantile(t_2,na.rm = T)
# max time in real data is 21.40 but we want a longer follow-up time
# we do not assume a distribution for censoring for simplification:
# no one leaves during the follow-up
T_max<-40
C<-T_max
# two clusters 
id_j.1<-sample(c(1,2),n_aid,replace = T,prob = c(0.65,0.35))
dat_event<-data.frame(cbind(id=1:n_aid,
                        population=id_j.1,t_1,t_2,C,x))
# combine three datasets
lmm_data<-cbind(y.pre,y.post[,c("y_1.post","y_2.post")])
lmm_data<-lmm_data%>%group_by(id)%>%mutate(population=id_j.1[id])%>%ungroup()
#remove the patients with "NA" in the event time
dat_event<-dat_event%>%na.omit()
lmm_data<-lmm_data%>%filter(id%in%dat_event$id)
#redefine id
dat_event$id<-1:nrow(dat_event)
lmm_data$id<-rep(c(1:nrow(dat_event)),each=5)
n_aid<-nrow(dat_event)
# select the observed time-to-event and y
dat_event<-dat_event%>%mutate(t=ifelse(population==1,t_1,t_2))%>%
  mutate(status=ifelse(t>C,0,1),
         event_time=ifelse(status==1,t,C))%>%select(-t)%>%
  dplyr::select(-t_1,-t_2,-C)
dat_combine.aid<-inner_join(lmm_data,dat_event,by="id")
dat_combine.aid<-dat_combine.aid%>%mutate(pop=population.x)%>%
  select(-population.x,-population.y)
dat_combine.aid<-dat_combine.aid%>%
  mutate(y1=ifelse(t<event_time,y_1.pre,y_1.post),
         y2=ifelse(t<event_time,y_2.pre,y_2.post))%>%
  select(-y_1.pre,-y_1.post,-y_2.pre,-y_2.post)
dat_combine.aid<-dat_combine.aid%>%
  mutate(y=ifelse(pop==1,y1,y2))%>%
  dplyr::select(-y1,-y2)
######no aids group##############
random_intercept1<-rnorm(n_no_aid,0,rand_int_sigma1)
random_intercept2<-rnorm(n_no_aid,0,rand_int_sigma2)
## pre-event
y.pre<-data.frame()
for (i in 1:n_no_aid) {
  y_1.pre<-t_no_missing*beta_t1.1+random_intercept1[i]+ #.1
    rnorm(length(t_no_missing),0,sigma_1)+Itcpt1
  y_2.pre<-t_no_missing*beta_t2.1+random_intercept2[i]+
    rnorm(length(t_no_missing),0,sigma_2)+Itcpt2
  y.pre<-rbind(y.pre,
               cbind(id=rep(i,length(t_no_missing)),y_1.pre,y_2.pre,
                     t=t_no_missing))
}
## post-event
y.post<-data.frame()
for (i in 1:n_no_aid) {
  y_1.post<-t_no_missing*beta_t1.2+random_intercept1[i]+ #.2
    rnorm(length(t_no_missing),0,sigma_1)+Itcpt1
  y_2.post<-t_no_missing*beta_t2.2+random_intercept2[i]+
    rnorm(length(t_no_missing),0,sigma_2)+Itcpt2
  y.post<-rbind(y.post,
                cbind(id=rep(i,length(t_no_missing)),y_1.post,y_2.post,
                      t=t_no_missing))
}
## event time
x<-sample(c(0,1),n_no_aid,replace = TRUE) #baseline covariate: drug
t_1<-c()
for (i in 1:n_no_aid) {
  t_1[i]<-newton_bisec(alpha=alpha1,X=x[i],b=random_intercept1[i],
                       beta=beta_d1,C=assoct1,a=beta_t1.1,intercept=Itcpt1,
                       eta=0,tol = 1e-10)      #.1
}
quantile(t_1,na.rm = T)
t_2<-c()
for (i in 1:n_no_aid) {
  t_2[i]<-newton_bisec(alpha=alpha2,X=x[i],b=random_intercept2[i],
                       beta=beta_d2,C=assoct2,a=beta_t2.1,intercept=Itcpt2,
                       eta=0,tol = 1e-10)      #.1
}
quantile(t_2,na.rm = T)
T_max<-40
C<-T_max
# two clusters 
id_j.2<-sample(c(1,2),n_no_aid,replace = T,prob = c(0.35,0.65))
dat_event<-data.frame(cbind(id=1:n_no_aid,
                            population=id_j.2,t_1,t_2,C,x))
# combine three datasets
lmm_data<-cbind(y.pre,y.post[,c("y_1.post","y_2.post")])
lmm_data<-lmm_data%>%group_by(id)%>%mutate(population=id_j.2[id])%>%ungroup()
#remove the patients with "NA" in the event time
dat_event<-dat_event%>%na.omit()
lmm_data<-lmm_data%>%filter(id%in%dat_event$id)
#redefine id
dat_event$id<-1:nrow(dat_event)
lmm_data$id<-rep(c(1:nrow(dat_event)),each=5)
n_no_aid<-nrow(dat_event)
# select the observed time-to-event and y
dat_event<-dat_event%>%mutate(t=ifelse(population==1,t_1,t_2))%>%
  mutate(status=ifelse(t>C,0,1),
         event_time=ifelse(status==1,t,C))%>%select(-t)%>%
dplyr::select(-t_1,-t_2,-C)
dat_combine.no.aid<-inner_join(lmm_data,dat_event,by="id")
dat_combine.no.aid<<-dat_combine.no.aid%>%mutate(pop=population.x)%>%
  select(-population.x,-population.y)
dat_combine.no.aid<-dat_combine.no.aid%>%
  mutate(y1=ifelse(t<event_time,y_1.pre,y_1.post),
         y2=ifelse(t<event_time,y_2.pre,y_2.post))%>%
  select(-y_1.pre,-y_1.post,-y_2.pre,-y_2.post)
dat_combine.no.aid<-dat_combine.no.aid%>%
  mutate(y=ifelse(pop==1,y1,y2))%>%
  dplyr::select(-y1,-y2)


######combine two sub-groups
dat<-rbind(cbind(data.frame(dat_combine.aid),J=1),
                   cbind(data.frame(dat_combine.no.aid),J=2))
dat<-dat%>%
  mutate(id=rep(1:(n_aid+n_no_aid),each=5)
         )
# range of CD4 in the real data is (0.00000 24.12468)
# 24 means 2400 cells/mm³
# normal range in real life 5-15 (500 - 1500)
range(dat$y)
# we remove the CD4 that is below than 0 which brings some missing data
dat<-dat%>%filter(y>=0)
# redefine id
id_left<-c(1:length(unique(dat$id)))
id_num<-dat%>%group_by(id)%>%summarize(n=n())
id_new<-c()
for (i in 1:length(id_left)) {
  id_new<-c(id_new,rep(id_left[i],id_num$n[i]))
}
dat$id<-id_new
save(dat,file = paste0("dat_test",sim.order,".rdata")) }


