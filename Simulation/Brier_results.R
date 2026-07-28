load("Brier_HDP_sim.rdata")
load("Brier_DP_Each_sim.rdata")
load("Brier_DP_All_sim.rdata")
load("Brier_Surv_sim.rdata")
load("Brier_JM_sim.rdata")
####t=2 u=3
brier.2.3<-c()
for (i in 1:100) {
  brier.2.3<-rbind(brier.2.3,
                   c(Brier_HDP_sim[[i]][1,1],Brier_DP_All_sim[[i]][1,1],
                     Brier_DP_Each_sim[[i]][1,1],Brier_Surv_sim[[i]][1,1],
                     Brier_JM_sim[[i]][1,1]))
}
brier.2.3<-data.frame(brier.2.3)
colnames(brier.2.3)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(brier.2.3)
sd(brier.2.3$HDP)
sd(brier.2.3$DPAll)
sd(brier.2.3$DPEach)
sd(brier.2.3$Surv)
sd(brier.2.3$JM)
wilcox.test(brier.2.3$HDP, brier.2.3$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.2.3$HDP, brier.2.3$DPEach,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.2.3$HDP, brier.2.3$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.2.3$HDP, brier.2.3$JM,
            paired = TRUE, exact = FALSE)#HDP best
####t=2 u=4
brier.2.4<-c()
for (i in 1:100) {
  brier.2.4<-rbind(brier.2.4,
                   c(Brier_HDP_sim[[i]][2,1],Brier_DP_All_sim[[i]][2,1],
                     Brier_DP_Each_sim[[i]][2,1],Brier_Surv_sim[[i]][2,1],
                     Brier_JM_sim[[i]][2,1]))
}
brier.2.4<-data.frame(brier.2.4)
colnames(brier.2.4)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(brier.2.4)
sd(brier.2.4$HDP)
sd(brier.2.4$DPAll)
sd(brier.2.4$DPEach)
sd(brier.2.4$Surv)
sd(brier.2.4$JM)
wilcox.test(brier.2.4$HDP, brier.2.4$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.2.4$HDP, brier.2.4$DPEach,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.2.4$HDP, brier.2.4$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.2.4$HDP, brier.2.4$JM,
            paired = TRUE, exact = FALSE) #HDP best
####t=2 u=5
brier.2.5<-c()
for (i in 1:100) {
  brier.2.5<-rbind(brier.2.5,
                   c(Brier_HDP_sim[[i]][3,1],Brier_DP_All_sim[[i]][3,1],
                     Brier_DP_Each_sim[[i]][3,1],Brier_Surv_sim[[i]][3,1],
                     Brier_JM_sim[[i]][3,1]))
}
brier.2.5<-data.frame(brier.2.5)
colnames(brier.2.5)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(brier.2.5)
sd(brier.2.5$HDP)
sd(brier.2.5$DPAll)
sd(brier.2.5$DPEach)
sd(brier.2.5$Surv)
sd(brier.2.5$JM)
wilcox.test(brier.2.5$HDP, brier.2.5$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.2.5$HDP, brier.2.5$DPEach,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.2.5$HDP, brier.2.5$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.2.5$HDP, brier.2.5$JM,
            paired = TRUE, exact = FALSE) #HDP best
####t=6 u=7
brier.6.7<-c()
for (i in 1:100) {
  brier.6.7<-rbind(brier.6.7,
                   c(Brier_HDP_sim[[i]][1,2],Brier_DP_All_sim[[i]][1,2],
                     Brier_DP_Each_sim[[i]][1,2],Brier_Surv_sim[[i]][1,2],
                     Brier_JM_sim[[i]][1,2]))
}
brier.6.7<-data.frame(brier.6.7)
colnames(brier.6.7)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(brier.6.7)
sd(brier.6.7$HDP)
sd(brier.6.7$DPAll)
sd(brier.6.7$DPEach)
sd(brier.6.7$Surv)
sd(brier.6.7$JM)
wilcox.test(brier.6.7$HDP, brier.6.7$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.6.7$HDP, brier.6.7$DPEach,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.6.7$HDP, brier.6.7$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.6.7$HDP, brier.6.7$JM,
            paired = TRUE, exact = FALSE)  #HDP best = DPAll
####t=6 u=8
brier.6.8<-c()
for (i in 1:100) {
  brier.6.8<-rbind(brier.6.8,
                   c(Brier_HDP_sim[[i]][2,2],Brier_DP_All_sim[[i]][2,2],
                     Brier_DP_Each_sim[[i]][2,2],Brier_Surv_sim[[i]][2,2],
                     Brier_JM_sim[[i]][2,2]))
}
brier.6.8<-data.frame(brier.6.8)
colnames(brier.6.8)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(brier.6.8)
sd(brier.6.8$HDP)
sd(brier.6.8$DPAll)
sd(brier.6.8$DPEach)
sd(brier.6.8$Surv)
sd(brier.6.8$JM)
wilcox.test(brier.6.8$HDP, brier.6.8$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.6.8$HDP, brier.6.8$DPEach,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.6.8$HDP, brier.6.8$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.6.8$HDP, brier.6.8$JM,
            paired = TRUE, exact = FALSE)#HDP best = DPAll

####t=6 u=9
brier.6.9<-c()
for (i in 1:100) {
  brier.6.9<-rbind(brier.6.9,
                   c(Brier_HDP_sim[[i]][3,2],Brier_DP_All_sim[[i]][3,2],
                     Brier_DP_Each_sim[[i]][3,2],Brier_Surv_sim[[i]][3,2],
                     Brier_JM_sim[[i]][3,2]))
}
brier.6.9<-data.frame(brier.6.9)
colnames(brier.6.9)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(brier.6.9)
sd(brier.6.9$HDP)
sd(brier.6.9$DPAll)
sd(brier.6.9$DPEach)
sd(brier.6.9$Surv)
sd(brier.6.9$JM)
wilcox.test(brier.6.9$HDP, brier.6.9$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.6.9$HDP, brier.6.9$DPEach,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.6.9$HDP, brier.6.9$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.6.9$HDP, brier.6.9$JM,
            paired = TRUE, exact = FALSE)
####t=12 u=13
brier.12.13<-c()
for (i in 1:100) {
  brier.12.13<-rbind(brier.12.13,
                   c(Brier_HDP_sim[[i]][1,3],Brier_DP_All_sim[[i]][1,3],
                     Brier_DP_Each_sim[[i]][1,3],Brier_Surv_sim[[i]][1,3],
                     Brier_JM_sim[[i]][1,3]))
}
brier.12.13<-data.frame(brier.12.13)
colnames(brier.12.13)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(brier.12.13)
sd(brier.12.13$HDP)
sd(brier.12.13$DPAll)
sd(brier.12.13$DPEach)
sd(brier.12.13$Surv)
sd(brier.12.13$JM)
wilcox.test(brier.12.13$HDP, brier.12.13$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.12.13$HDP, brier.12.13$DPEach,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.12.13$HDP, brier.12.13$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.12.13$HDP, brier.12.13$JM,
            paired = TRUE, exact = FALSE)#HDP best = DPAll
####t=12 u=14
brier.12.14<-c()
for (i in 1:100) {
  brier.12.14<-rbind(brier.12.14,
                   c(Brier_HDP_sim[[i]][2,3],Brier_DP_All_sim[[i]][2,3],
                     Brier_DP_Each_sim[[i]][2,3],Brier_Surv_sim[[i]][2,3],
                     Brier_JM_sim[[i]][2,3]))
}
brier.12.14<-data.frame(brier.12.14)
colnames(brier.12.14)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(brier.12.14)
sd(brier.12.14$HDP)
sd(brier.12.14$DPAll)
sd(brier.12.14$DPEach)
sd(brier.12.14$Surv)
sd(brier.12.14$JM)
wilcox.test(brier.12.14$HDP, brier.12.14$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.12.14$HDP, brier.12.14$DPEach,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.12.14$HDP, brier.12.14$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.12.14$HDP, brier.12.14$JM,
            paired = TRUE, exact = FALSE)#HDP best

####t=12 u=15
brier.12.15<-c()
for (i in 1:100) {
  brier.12.15<-rbind(brier.12.15,
                   c(Brier_HDP_sim[[i]][3,3],Brier_DP_All_sim[[i]][3,3],
                     Brier_DP_Each_sim[[i]][3,3],Brier_Surv_sim[[i]][3,3],
                     Brier_JM_sim[[i]][3,3]))
}
brier.12.15<-data.frame(brier.12.15)
colnames(brier.12.15)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(brier.12.15)
sd(brier.12.15$HDP)
sd(brier.12.15$DPAll)
sd(brier.12.15$DPEach)
sd(brier.12.15$Surv)
sd(brier.12.15$JM)
wilcox.test(brier.12.15$HDP, brier.12.15$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.12.15$HDP, brier.12.15$DPEach,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.12.15$HDP, brier.12.15$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.12.15$HDP, brier.12.15$JM,
            paired = TRUE, exact = FALSE) #HDP best = DPAll

####t=18 u=19
brier.18.19<-c()
for (i in 1:100) {
  brier.18.19<-rbind(brier.18.19,
                     c(Brier_HDP_sim[[i]][1,4],Brier_DP_All_sim[[i]][1,4],
                       Brier_DP_Each_sim[[i]][1,4],Brier_Surv_sim[[i]][1,4],
                       Brier_JM_sim[[i]][1,4]))
}
brier.18.19<-data.frame(brier.18.19)
colnames(brier.18.19)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(brier.18.19)
sd(brier.18.19$HDP)
sd(brier.18.19$DPAll)
sd(brier.18.19$DPEach)
sd(brier.18.19$Surv)
sd(brier.18.19$JM)
wilcox.test(brier.18.19$HDP, brier.18.19$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.18.19$HDP, brier.18.19$DPEach,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.18.19$HDP, brier.18.19$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.18.19$HDP, brier.18.19$JM,
            paired = TRUE, exact = FALSE)  #HDP best = DPAll

####t=18 u=20
brier.18.20<-c()
for (i in 1:100) {
  brier.18.20<-rbind(brier.18.20,
                     c(Brier_HDP_sim[[i]][2,4],Brier_DP_All_sim[[i]][2,4],
                       Brier_DP_Each_sim[[i]][2,4],Brier_Surv_sim[[i]][2,4],
                       Brier_JM_sim[[i]][2,4]))
}
brier.18.20<-data.frame(brier.18.20)
colnames(brier.18.20)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(brier.18.20)
sd(brier.18.20$HDP)
sd(brier.18.20$DPAll)
sd(brier.18.20$DPEach)
sd(brier.18.20$Surv)
sd(brier.18.20$JM)
wilcox.test(brier.18.20$HDP, brier.18.20$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.18.20$HDP, brier.18.20$DPEach,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.18.20$HDP, brier.18.20$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.18.20$HDP, brier.18.20$JM,
            paired = TRUE, exact = FALSE)   #HDP best = DPAll

####t=18 u=21
brier.18.21<-c()
for (i in 1:100) {
  brier.18.21<-rbind(brier.18.21,
                     c(Brier_HDP_sim[[i]][3,4],Brier_DP_All_sim[[i]][3,4],
                       Brier_DP_Each_sim[[i]][3,4],Brier_Surv_sim[[i]][3,4],
                       Brier_JM_sim[[i]][3,4]))
}
brier.18.21<-data.frame(brier.18.21)
colnames(brier.18.21)<-c("HDP","DPAll","DPEach","Surv","JM")
colMeans(brier.18.21)
sd(brier.18.21$HDP)
sd(brier.18.21$DPAll)
sd(brier.18.21$DPEach)
sd(brier.18.21$Surv)
sd(brier.18.21$JM)
wilcox.test(brier.18.21$HDP, brier.18.21$DPAll,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.18.21$HDP, brier.18.21$DPEach,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.18.21$HDP, brier.18.21$Surv,
            paired = TRUE, exact = FALSE)
wilcox.test(brier.18.21$HDP, brier.18.21$JM,
            paired = TRUE, exact = FALSE)
