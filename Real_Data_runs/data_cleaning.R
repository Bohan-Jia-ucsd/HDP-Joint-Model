library(tidyverse)
library(lubridate)
library(ggplot2)
library(R2jags)
# 594
cd4<-read.csv("CD4Log.csv")
###check missing
sum(is.na(cd4$cd4per))
sum(is.na(cd4$cd4abs)) #cell/mm^3
sum(is.na(cd4$cd8per))
sum(is.na(cd4$cd8))
longitudinal594<-cd4%>%
  select(pid, visit_date, cd4abs)
longitudinal594<-longitudinal594[order(longitudinal594$pid,longitudinal594$visit_date),]
longitudinal594<-longitudinal594%>%na.omit()
longitudinal594$visit_date<-as.Date(longitudinal594$visit_date)
# from paper
start_date<-as.Date("2014-01-01")
start_date<-as.numeric(start_date)
longitudinal594$visit_date<-as.numeric(longitudinal594$visit_date)
longitudinal594$visit_date<-longitudinal594$visit_date-start_date
ggplot()+
  geom_point(data = longitudinal594,aes(x=visit_date,y=cd4abs,group=pid,color=pid))+
  geom_line(data = longitudinal594,aes(x=visit_date,y=cd4abs,group=pid,color=pid))+
  theme(legend.position = "none")+
  geom_smooth(data = longitudinal594,aes(x=visit_date,y=cd4abs,group=FALSE),
              color="black")
unique(longitudinal594$pid)
#will add more baseline covariates into this dataset

## time-to-event
arv594<-read.csv("ARVLog.csv")
sum(arv594$pid%in%unique(longitudinal594$pid))
arv594[arv594$start_date=="","start_date"]<-"2016-03-31"
arv594$start_date<-as.numeric(as.Date(arv594$start_date))-start_date
arv594<-arv594[order(arv594$pid,arv594$start_date),]
table(arv594$start_date)
arv5942<-arv594%>%select(pid,start_date)%>%
  group_by(pid)%>%slice(1)%>%ungroup()
table(arv5942$start_date)
arv5942<-arv5942%>%filter(start_date>0)
sum(arv5942$pid%in%unique(longitudinal594$pid))

dat_combine_complete594<-inner_join(longitudinal594,arv5942,by="pid")

dat_combine_complete594<-dat_combine_complete594%>%
  mutate(status=ifelse(start_date==820,0,1))
dat_surv594<-dat_combine_complete594%>%group_by(pid)%>%slice(1)%>%ungroup()%>%
  select(pid,start_date,status)
dat_surv594<-data.frame(dat_surv594)
range(dat_surv594$start_date)
range(dat_combine_complete594$visit_date)
library(survival)
library(survminer)
fit_surv<-survfit(Surv(start_date, status) ~ 1, data = dat_surv594)
ggsurvplot(fit_surv, data = dat_surv594)
### re-scale time
# have to make sure the event times and the CD4 times are on the same scale and
# have to make sure two dataset have the same time scale
range(dat_combine_complete594$visit_date)
range(dat_surv594$start_date)
dat_surv594$start_date<-dat_surv594$start_date/820
dat_combine_complete594$visit_date<-dat_combine_complete594$visit_date/820
fit_surv<-survfit(Surv(start_date, status) ~ 1, data = dat_surv594)
ggsurvplot(fit_surv, data = dat_surv594)
ggplot()+
  geom_point(data = dat_combine_complete594,aes(x=visit_date,y=cd4abs,group=pid,color=pid))+
  geom_line(data = dat_combine_complete594,aes(x=visit_date,y=cd4abs,group=pid,color=pid))+
  theme(legend.position = "none")+
  geom_smooth(data = dat_combine_complete594,aes(x=visit_date,y=cd4abs,group=FALSE),
              color="black")
# 320
datCD4 <- read.table("ACTG320_CD4.txt", header = TRUE) #cell/mm^3 from paper
datARV <- read.table("ACTG320_RX.txt", header = TRUE)
range(datCD4$CD4Date)
datCD4<-datCD4%>%filter(CD4Date>=0)
datARV<-datARV%>%filter(StartDate>0)
datARV<-datARV%>%group_by(PtID)%>%slice(1)%>%ungroup()
ID_left<-datARV$PtID
datCD4<-datCD4%>%filter(PtID %in% ID_left)
sum(unique(datCD4$PtID)==datARV$PtID)
datCD4$PtID<-as.factor(datCD4$PtID)
ggplot()+
  geom_point(data = datCD4,aes(x=CD4Date,y=CD4Count,group=PtID,color=PtID))+
  geom_line(data = datCD4,aes(x=CD4Date,y=CD4Count,group=PtID,color=PtID))+
  theme(legend.position = "none")+
  geom_smooth(data = datCD4,aes(x=CD4Date,y=CD4Count,group=FALSE),
              color="black")
datARV<-datARV%>%mutate(status=1)
fit_surv<-survfit(Surv(StartDate, status) ~ 1, data = datARV)
ggsurvplot(fit_surv, data = datARV)
dat_surv320<-datARV%>%mutate(start_date=StartDate)%>%
  select(PtID,start_date,status)
dat_combine_complete320<-datCD4
#re-scale time
range(dat_combine_complete320$CD4Date)
range(dat_surv320$start_date)
dat_surv320$start_date<-dat_surv320$start_date/56
dat_combine_complete320$CD4Date<-dat_combine_complete320$CD4Date/56
ggplot()+
  geom_point(data = dat_combine_complete320,aes(x=CD4Date,y=CD4Count,group=PtID,color=PtID))+
  geom_line(data = dat_combine_complete320,aes(x=CD4Date,y=CD4Count,group=PtID,color=PtID))+
  theme(legend.position = "none")+
  geom_smooth(data = dat_combine_complete320,aes(x=CD4Date,y=CD4Count,group=FALSE),
              color="black")  
fit_surv<-survfit(Surv(start_date, status) ~ 1, data = dat_surv320)
ggsurvplot(fit_surv, data = dat_surv320)
# CD4 cell
range(dat_combine_complete594$cd4abs) #cell/mm^3
range(dat_combine_complete320$CD4Count) #cell/mm^3 from paper
# re-scale them
dat_combine_complete594$cd4abs<-dat_combine_complete594$cd4abs/100
dat_combine_complete320$CD4Count<-dat_combine_complete320$CD4Count/100

#########combine two population

dat_combine_complete320<-dat_combine_complete320%>%mutate(id=PtID,visit=CD4Date,
                                                          cd4count=CD4Count,
                                                          pop="320")%>%
  select(-PtID,-CD4Date,-CD4Count)
dat_combine_complete594<-dat_combine_complete594%>%mutate(id=pid,visit=visit_date,
                                                          cd4count=cd4abs,
                                                          pop="594")%>%
  select(id,visit,cd4count,pop)
dat_combine_complete320<-data.frame(dat_combine_complete320)
dat_combine_complete594<-data.frame(dat_combine_complete594)
dat_CD4_combine<-rbind(dat_combine_complete594,dat_combine_complete320)  

dat_surv320<-data.frame(dat_surv320)
dat_surv320<-dat_surv320%>%mutate(id=PtID,ART_init=start_date,
                                  pop="320")%>%select(id,ART_init,status,pop)
dat_surv594<-data.frame(dat_surv594)
dat_surv594<-dat_surv594%>%mutate(id=pid,ART_init=start_date,
                                  pop="594")%>%select(id,ART_init,status,pop)
dat_surv_combine<-rbind(dat_surv594,dat_surv320) 

ggplot()+
  geom_point(data = dat_CD4_combine,aes(x=visit,y=cd4count,group=id,color=pop))+
  geom_line(data = dat_CD4_combine,aes(x=visit,y=cd4count,group=id,color=pop))+
  ylab("CD4 Count")+xlab("Time")

fit_surv<-survfit(Surv(ART_init, status) ~ pop, data = dat_surv_combine)
ggsurvplot(fit_surv, data = dat_surv_combine)


km_df <- broom::tidy(fit_surv) %>%
  filter(time > 0, estimate > 0, estimate < 1) %>%
  mutate(
    log_time = log(time),
    log_neglog_surv = log(-log(estimate))
  )

ggplot(km_df, aes(x = log_time, y = log_neglog_surv, color = strata)) +
  geom_step(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "log(Time)",
    y = "log[-log(KM survival)]",
    color = "Cohort"
  ) +
  theme_bw()

dat_surv_combine$id==unique(dat_CD4_combine$id)
id_num<-dat_CD4_combine%>%group_by(id)%>%summarize(num=n())
dat_surv_combine$id<-1:nrow(dat_surv_combine)
id_new<-c()
for (i in 1:nrow(dat_surv_combine)) {
  id_new<-c(id_new,rep(i,id_num$num[i]))
}
dat_CD4_combine$id<-id_new
save(dat_CD4_combine, file = paste0("dat_CD4_combine.rdata"))
save(dat_surv_combine, file = paste0("dat_surv_combine.rdata"))
load(paste0("dat_CD4_combine.rdata"))
load(paste0("dat_surv_combine.rdata"))
