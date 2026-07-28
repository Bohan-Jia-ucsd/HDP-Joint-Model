load("AUC_HDP_sim.rdata")
load("AUC_DP_Each_sim.rdata")
load("AUC_DP_All_sim.rdata")
load("AUC_Surv_sim.rdata")
load("AUC_JM_sim.rdata")
library(pracma)
library(tidyverse)
####t=2 u=3
AUC.2.3<-c()
for (i in 1:100) {
  AUC.2.3<-rbind(AUC.2.3,
                   c(trapz(1-AUC_HDP_sim[[i]][,3], AUC_HDP_sim[[i]][,2]),
                     trapz(1-AUC_DP_All_sim[[i]][,3], AUC_DP_All_sim[[i]][,2]),
                     trapz(1-AUC_DP_Each_sim[[i]][,3], AUC_DP_Each_sim[[i]][,2]),
                     trapz(1-AUC_Surv_sim[[i]][,3], AUC_Surv_sim[[i]][,2]),
                     trapz(1-AUC_JM_sim[[i]][,3], AUC_JM_sim[[i]][,2])))
}
AUC.2.3<-data.frame(AUC.2.3)
AUC.2.3<-na.omit(AUC.2.3) 
# Read Please
# reason why we have "na" here is there is no event time in test data73
# which is people who survives at time 2 also survives at time 3
# pretty normal because 2 and 3 are too close to each other!
colnames(AUC.2.3)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(AUC.2.3)
sd(AUC.2.3$HDP)
sd(AUC.2.3$DPAll)
sd(AUC.2.3$DPEach)
sd(AUC.2.3$Surv)
sd(AUC.2.3$JM)
wilcox.test(AUC.2.3$HDP, AUC.2.3$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.2.3$HDP, AUC.2.3$DPEach,
            paired = TRUE, exact = FALSE) 
wilcox.test(AUC.2.3$HDP, AUC.2.3$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.2.3$HDP, AUC.2.3$JM,
            paired = TRUE, exact = FALSE)#HDP best
####t=2 u=4
AUC.2.4<-c()
for (i in 1:100) {
  AUC.2.4<-rbind(AUC.2.4,
                 c(trapz(1-AUC_HDP_sim[[i]][,5], AUC_HDP_sim[[i]][,4]),
                   trapz(1-AUC_DP_All_sim[[i]][,5], AUC_DP_All_sim[[i]][,4]),
                   trapz(1-AUC_DP_Each_sim[[i]][,5], AUC_DP_Each_sim[[i]][,4]),
                   trapz(1-AUC_Surv_sim[[i]][,5], AUC_Surv_sim[[i]][,4]),
                   trapz(1-AUC_JM_sim[[i]][,5], AUC_JM_sim[[i]][,4])))
}
AUC.2.4<-data.frame(AUC.2.4)
colnames(AUC.2.4)<-c("HDP","DPAll","DPEach","Surv","JM")
AUC.2.4<-na.omit(AUC.2.4) 
colMeans(AUC.2.4)
sd(AUC.2.4$HDP)
sd(AUC.2.4$DPAll)
sd(AUC.2.4$DPEach)
sd(AUC.2.4$Surv)
sd(AUC.2.4$JM)
wilcox.test(AUC.2.4$HDP, AUC.2.4$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.2.4$HDP, AUC.2.4$DPEach,
            paired = TRUE, exact = FALSE) 
wilcox.test(AUC.2.4$HDP, AUC.2.4$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.2.4$HDP, AUC.2.4$JM,
            paired = TRUE, exact = FALSE)#HDP best
#HDP best
####t=2 u=5
AUC.2.5<-c()
for (i in 1:100) {
  AUC.2.5<-rbind(AUC.2.5,
                 c(trapz(1-AUC_HDP_sim[[i]][,7], AUC_HDP_sim[[i]][,6]),
                   trapz(1-AUC_DP_All_sim[[i]][,7], AUC_DP_All_sim[[i]][,6]),
                   trapz(1-AUC_DP_Each_sim[[i]][,7], AUC_DP_Each_sim[[i]][,6]),
                   trapz(1-AUC_Surv_sim[[i]][,7], AUC_Surv_sim[[i]][,6]),
                   trapz(1-AUC_JM_sim[[i]][,7], AUC_JM_sim[[i]][,6])))
}
AUC.2.5<-data.frame(AUC.2.5)
colnames(AUC.2.5)<-c("HDP","DPAll","DPEach","Surv","JM")
AUC.2.5<-na.omit(AUC.2.5) 
colMeans(AUC.2.5)
sd(AUC.2.5$HDP)
sd(AUC.2.5$DPAll)
sd(AUC.2.5$DPEach)
sd(AUC.2.5$Surv)
sd(AUC.2.5$JM)
wilcox.test(AUC.2.5$HDP, AUC.2.5$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.2.5$HDP, AUC.2.5$DPEach,
            paired = TRUE, exact = FALSE) 
wilcox.test(AUC.2.5$HDP, AUC.2.5$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.2.5$HDP, AUC.2.5$JM,
            paired = TRUE, exact = FALSE) #HDP best
####t=6 u=7
AUC.6.7<-c()
for (i in 1:100) {
  AUC.6.7<-rbind(AUC.6.7,
                 c(trapz(1-AUC_HDP_sim[[i]][,9], AUC_HDP_sim[[i]][,8]),
                   trapz(1-AUC_DP_All_sim[[i]][,9], AUC_DP_All_sim[[i]][,8]),
                   trapz(1-AUC_DP_Each_sim[[i]][,9], AUC_DP_Each_sim[[i]][,8]),
                   trapz(1-AUC_Surv_sim[[i]][,9], AUC_Surv_sim[[i]][,8]),
                   trapz(1-AUC_JM_sim[[i]][,9], AUC_JM_sim[[i]][,8])))
}
AUC.6.7<-data.frame(AUC.6.7)
colnames(AUC.6.7)<-c("HDP","DPAll","DPEach","Surv","JM")
AUC.6.7<-na.omit(AUC.6.7) 
colMeans(AUC.6.7)
sd(AUC.6.7$HDP)
sd(AUC.6.7$DPAll)
sd(AUC.6.7$DPEach)
sd(AUC.6.7$Surv)
sd(AUC.6.7$JM)
wilcox.test(AUC.6.7$HDP, AUC.6.7$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.6.7$HDP, AUC.6.7$DPEach,
            paired = TRUE, exact = FALSE) 
wilcox.test(AUC.6.7$HDP, AUC.6.7$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.6.7$HDP, AUC.6.7$JM,
            paired = TRUE, exact = FALSE) #HDP best = DPAll
####t=6 u=8
AUC.6.8<-c()
for (i in 1:100) {
  AUC.6.8<-rbind(AUC.6.8,
                 c(trapz(1-AUC_HDP_sim[[i]][,11], AUC_HDP_sim[[i]][,10]),
                   trapz(1-AUC_DP_All_sim[[i]][,11], AUC_DP_All_sim[[i]][,10]),
                   trapz(1-AUC_DP_Each_sim[[i]][,11], AUC_DP_Each_sim[[i]][,10]),
                   trapz(1-AUC_Surv_sim[[i]][,11], AUC_Surv_sim[[i]][,10]),
                   trapz(1-AUC_JM_sim[[i]][,11], AUC_JM_sim[[i]][,10])))
}
AUC.6.8<-data.frame(AUC.6.8)
colnames(AUC.6.8)<-c("HDP","DPAll","DPEach","Surv","JM")
AUC.6.8<-na.omit(AUC.6.8) 
colMeans(AUC.6.8)
sd(AUC.6.8$HDP)
sd(AUC.6.8$DPAll)
sd(AUC.6.8$DPEach)
sd(AUC.6.8$Surv)
sd(AUC.6.8$JM)
wilcox.test(AUC.6.8$HDP, AUC.6.8$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.6.8$HDP, AUC.6.8$DPEach,
            paired = TRUE, exact = FALSE) 
wilcox.test(AUC.6.8$HDP, AUC.6.8$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.6.8$HDP, AUC.6.8$JM,
            paired = TRUE, exact = FALSE)#HDP best = DPAll

####t=6 u=9
AUC.6.9<-c()
for (i in 1:100) {
  AUC.6.9<-rbind(AUC.6.9,
                 c(trapz(1-AUC_HDP_sim[[i]][,13], AUC_HDP_sim[[i]][,12]),
                   trapz(1-AUC_DP_All_sim[[i]][,13], AUC_DP_All_sim[[i]][,12]),
                   trapz(1-AUC_DP_Each_sim[[i]][,13], AUC_DP_Each_sim[[i]][,12]),
                   trapz(1-AUC_Surv_sim[[i]][,13], AUC_Surv_sim[[i]][,12]),
                   trapz(1-AUC_JM_sim[[i]][,13], AUC_JM_sim[[i]][,12])))
}
AUC.6.9<-data.frame(AUC.6.9)
colnames(AUC.6.9)<-c("HDP","DPAll","DPEach","Surv","JM")
AUC.6.9<-na.omit(AUC.6.9)
colMeans(AUC.6.9)
sd(AUC.6.9$HDP)
sd(AUC.6.9$DPAll)
sd(AUC.6.9$DPEach)
sd(AUC.6.9$Surv)
sd(AUC.6.9$JM)
wilcox.test(AUC.6.9$HDP, AUC.6.9$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.6.9$HDP, AUC.6.9$DPEach,
            paired = TRUE, exact = FALSE) 
wilcox.test(AUC.6.9$HDP, AUC.6.9$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.6.9$HDP, AUC.6.9$JM,
            paired = TRUE, exact = FALSE)#HDP best = DPAll

####t=12 u=13
AUC.12.13<-c()
for (i in 1:100) {
  AUC.12.13<-rbind(AUC.12.13,
                   c(trapz(1-AUC_HDP_sim[[i]][,15], AUC_HDP_sim[[i]][,14]),
                     trapz(1-AUC_DP_All_sim[[i]][,15], AUC_DP_All_sim[[i]][,14]),
                     trapz(1-AUC_DP_Each_sim[[i]][,15], AUC_DP_Each_sim[[i]][,14]),
                     trapz(1-AUC_Surv_sim[[i]][,15], AUC_Surv_sim[[i]][,14]),
                     trapz(1-AUC_JM_sim[[i]][,15], AUC_JM_sim[[i]][,14])))
}
AUC.12.13<-data.frame(AUC.12.13)
colnames(AUC.12.13)<-c("HDP","DPAll","DPEach","Surv","JM")
AUC.12.13<-na.omit(AUC.12.13)
colMeans(AUC.12.13)
sd(AUC.12.13$HDP)
sd(AUC.12.13$DPAll)
sd(AUC.12.13$DPEach)
sd(AUC.12.13$Surv)
sd(AUC.12.13$JM)
wilcox.test(AUC.12.13$HDP, AUC.12.13$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.12.13$HDP, AUC.12.13$DPEach,
            paired = TRUE, exact = FALSE) 
wilcox.test(AUC.12.13$HDP, AUC.12.13$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.12.13$HDP, AUC.12.13$JM,
            paired = TRUE, exact = FALSE) #HDP best = DPAll
####t=12 u=14
AUC.12.14<-c()
for (i in 1:100) {
  AUC.12.14<-rbind(AUC.12.14,
                   c(trapz(1-AUC_HDP_sim[[i]][,17], AUC_HDP_sim[[i]][,16]),
                     trapz(1-AUC_DP_All_sim[[i]][,17], AUC_DP_All_sim[[i]][,16]),
                     trapz(1-AUC_DP_Each_sim[[i]][,17], AUC_DP_Each_sim[[i]][,16]),
                     trapz(1-AUC_Surv_sim[[i]][,17], AUC_Surv_sim[[i]][,16]),
                     trapz(1-AUC_JM_sim[[i]][,17], AUC_JM_sim[[i]][,16])))
}
AUC.12.14<-data.frame(AUC.12.14)
colnames(AUC.12.14)<-c("HDP","DPAll","DPEach","Surv","JM")
AUC.12.14<-na.omit(AUC.12.14)
colMeans(AUC.12.14)
sd(AUC.12.14$HDP)
sd(AUC.12.14$DPAll)
sd(AUC.12.14$DPEach)
sd(AUC.12.14$Surv)
sd(AUC.12.14$JM)
wilcox.test(AUC.12.14$HDP, AUC.12.14$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.12.14$HDP, AUC.12.14$DPEach,
            paired = TRUE, exact = FALSE) 
wilcox.test(AUC.12.14$HDP, AUC.12.14$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.12.14$HDP, AUC.12.14$JM,
            paired = TRUE, exact = FALSE)  #HDP best

####t=12 u=15
AUC.12.15<-c()
for (i in 1:100) {
  AUC.12.15<-rbind(AUC.12.15,
                   c(trapz(1-AUC_HDP_sim[[i]][,19], AUC_HDP_sim[[i]][,18]),
                     trapz(1-AUC_DP_All_sim[[i]][,19], AUC_DP_All_sim[[i]][,18]),
                     trapz(1-AUC_DP_Each_sim[[i]][,19], AUC_DP_Each_sim[[i]][,18]),
                     trapz(1-AUC_Surv_sim[[i]][,19], AUC_Surv_sim[[i]][,18]),
                     trapz(1-AUC_JM_sim[[i]][,19], AUC_JM_sim[[i]][,18])))
}
AUC.12.15<-data.frame(AUC.12.15)
colnames(AUC.12.15)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(AUC.12.15)
sd(AUC.12.15$HDP)
sd(AUC.12.15$DPAll)
sd(AUC.12.15$DPEach)
sd(AUC.12.15$Surv)
sd(AUC.12.15$JM)
wilcox.test(AUC.12.15$HDP, AUC.12.15$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.12.15$HDP, AUC.12.15$DPEach,
            paired = TRUE, exact = FALSE) 
wilcox.test(AUC.12.15$HDP, AUC.12.15$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.12.15$HDP, AUC.12.15$JM,
            paired = TRUE, exact = FALSE)#HDP best

####t=18 u=19
AUC.18.19<-c()
for (i in 1:100) {
  AUC.18.19<-rbind(AUC.18.19,
                   c(trapz(1-AUC_HDP_sim[[i]][,21], AUC_HDP_sim[[i]][,20]),
                     trapz(1-AUC_DP_All_sim[[i]][,21], AUC_DP_All_sim[[i]][,20]),
                     trapz(1-AUC_DP_Each_sim[[i]][,21], AUC_DP_Each_sim[[i]][,20]),
                     trapz(1-AUC_Surv_sim[[i]][,21], AUC_Surv_sim[[i]][,20]),
                     trapz(1-AUC_JM_sim[[i]][,21], AUC_JM_sim[[i]][,20])))
}
AUC.18.19<-data.frame(AUC.18.19)
colnames(AUC.18.19)<-c("HDP","DPAll","DPEach","Surv","JM")
AUC.18.19<-na.omit(AUC.18.19)
colMeans(AUC.18.19)
sd(AUC.18.19$HDP)
sd(AUC.18.19$DPAll)
sd(AUC.18.19$DPEach)
sd(AUC.18.19$Surv)
sd(AUC.18.19$JM)
wilcox.test(AUC.18.19$HDP, AUC.18.19$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.18.19$HDP, AUC.18.19$DPEach,
            paired = TRUE, exact = FALSE) 
wilcox.test(AUC.18.19$HDP, AUC.18.19$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.18.19$HDP, AUC.18.19$JM,
            paired = TRUE, exact = FALSE) #HDP best = DPAll

####t=18 u=20
AUC.18.20<-c()
for (i in 1:100) {
  AUC.18.20<-rbind(AUC.18.20,
                   c(trapz(1-AUC_HDP_sim[[i]][,23], AUC_HDP_sim[[i]][,22]),
                     trapz(1-AUC_DP_All_sim[[i]][,23], AUC_DP_All_sim[[i]][,22]),
                     trapz(1-AUC_DP_Each_sim[[i]][,23], AUC_DP_Each_sim[[i]][,22]),
                     trapz(1-AUC_Surv_sim[[i]][,23], AUC_Surv_sim[[i]][,22]),
                     trapz(1-AUC_JM_sim[[i]][,23], AUC_JM_sim[[i]][,22])))
}
AUC.18.20<-data.frame(AUC.18.20)
colnames(AUC.18.20)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(AUC.18.20)
sd(AUC.18.20$HDP)
sd(AUC.18.20$DPAll)
sd(AUC.18.20$DPEach)
sd(AUC.18.20$Surv)
sd(AUC.18.20$JM)
wilcox.test(AUC.18.20$HDP, AUC.18.20$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.18.20$HDP, AUC.18.20$DPEach,
            paired = TRUE, exact = FALSE) 
wilcox.test(AUC.18.20$HDP, AUC.18.20$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.18.20$HDP, AUC.18.20$JM,
            paired = TRUE, exact = FALSE) #HDP best 

####t=18 u=21
AUC.18.21<-c()
for (i in 1:100) {
  AUC.18.21<-rbind(AUC.18.21,
                   c(trapz(1-AUC_HDP_sim[[i]][,25], AUC_HDP_sim[[i]][,24]),
                     trapz(1-AUC_DP_All_sim[[i]][,25], AUC_DP_All_sim[[i]][,24]),
                     trapz(1-AUC_DP_Each_sim[[i]][,25], AUC_DP_Each_sim[[i]][,24]),
                     trapz(1-AUC_Surv_sim[[i]][,25], AUC_Surv_sim[[i]][,24]),
                     trapz(1-AUC_JM_sim[[i]][,25], AUC_JM_sim[[i]][,24])))
}
AUC.18.21<-data.frame(AUC.18.21)
colnames(AUC.18.21)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(AUC.18.21)
sd(AUC.18.21$HDP)
sd(AUC.18.21$DPAll)
sd(AUC.18.21$DPEach)
sd(AUC.18.21$Surv)
sd(AUC.18.21$JM)
wilcox.test(AUC.18.21$HDP, AUC.18.21$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.18.21$HDP, AUC.18.21$DPEach,
            paired = TRUE, exact = FALSE) 
wilcox.test(AUC.18.21$HDP, AUC.18.21$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(AUC.18.21$HDP, AUC.18.21$JM,
            paired = TRUE, exact = FALSE)#HDP best
