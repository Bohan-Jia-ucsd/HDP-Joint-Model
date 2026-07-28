library(tidyverse)
library(fclust)
ARI_fuzzy.all<-rep(0,100)
ARI_fuzzy.1<-rep(0,100)
ARI_fuzzy.2<-rep(0,100)
ARI_Hard.all<-rep(0,100)
ARI_Hard.1<-rep(0,100)
ARI_Hard.2<-rep(0,100)
#HDP_joint
for (i in 1:100) {
load(paste0("Hier_weibull_sim_Comp",i,".rdata"))

load(paste0("sim_100_train/dat_train",i,".rdata"))
dat_lmm<-dat%>%select(-event_time,-status)
dat_surv<-dat%>%group_by(id)%>%slice(1)%>%ungroup()%>%
  select(-y,-t)
zeta_HDP<-as.data.frame(t(jags.m$BUGSoutput$sims.list$zeta[1:2000,]))
df_long <- zeta_HDP %>%
  mutate(item = row_number()) %>%
  pivot_longer(-item, names_to = "sample", values_to = "class")
prob_df <- df_long %>%
  group_by(item, class) %>%
  summarise(prob = n() / 2000, .groups = "drop")
prob_matrix <- prob_df %>%
  pivot_wider(names_from = class, values_from = prob, values_fill = 0)
#fuzzy
ARI_fuzzy.all[i]<-ARI.F(dat_surv$pop,prob_matrix[,-1])
ARI_fuzzy.1[i]<-ARI.F(dat_surv$pop[dat_surv$J==1],prob_matrix[dat_surv$J==1,-1])
ARI_fuzzy.2[i]<-ARI.F(dat_surv$pop[dat_surv$J==2],prob_matrix[dat_surv$J==2,-1])
#hard
prob_df <- zeta_HDP %>%
  mutate(item = row_number()) %>%
  pivot_longer(-item, values_to = "class") %>%
  count(item, class) %>%
  group_by(item) %>%
  mutate(prob = n / sum(n)) %>%
  select(-n)
hard_prob <- prob_df %>%
  group_by(item) %>%
  mutate(prob = if_else(prob == max(prob), 1, 0)) %>%
  ungroup()
hard_matrix <- hard_prob %>%
  pivot_wider(names_from = class, values_from = prob, values_fill = 0)
ARI_Hard.all[i]<-ARI.F(dat_surv$pop,hard_matrix[,-1])
ARI_Hard.1[i]<-ARI.F(dat_surv$pop[dat_surv$J==1],hard_matrix[dat_surv$J==1,-1])
ARI_Hard.2[i]<-ARI.F(dat_surv$pop[dat_surv$J==2],hard_matrix[dat_surv$J==2,-1])
}
ARI_HDP<-data.frame(cbind(ARI_fuzzy.all,ARI_fuzzy.1,ARI_fuzzy.2,
               ARI_Hard.all,ARI_Hard.1,ARI_Hard.2))

###################DP_ALL
sim_seed<-c(1:100)
ARI_fuzzy.all<-rep(0,length(sim_seed))
ARI_fuzzy.1<-rep(0,length(sim_seed))
ARI_fuzzy.2<-rep(0,length(sim_seed))
ARI_Hard.all<-rep(0,length(sim_seed))
ARI_Hard.1<-rep(0,length(sim_seed))
ARI_Hard.2<-rep(0,length(sim_seed))
for (i in 1:length(sim_seed)) {
load(paste0("DP_all_sim_Comp",sim_seed[i],".rdata"))
load(paste0("/Users/bohanjia/Desktop/research/gsr/Sim_Final/sim_100_train/dat_train",sim_seed[i],".rdata"))
dat_lmm<-dat%>%select(-event_time,-status)
dat_surv<-dat%>%group_by(id)%>%slice(1)%>%ungroup()%>%
    select(-y,-t)
zeta_DP<-as.data.frame(t(jags.m$BUGSoutput$sims.list$zeta[1:2000,]))
df_long <- zeta_DP %>%
  mutate(item = row_number()) %>%
  pivot_longer(-item, names_to = "sample", values_to = "class")
prob_df <- df_long %>%
  group_by(item, class) %>%
  summarise(prob = n() / 2000, .groups = "drop")
prob_matrix <- prob_df %>%
  pivot_wider(names_from = class, values_from = prob, values_fill = 0)
#fuzzy
ARI_fuzzy.all[i]<-ARI.F(dat_surv$pop,prob_matrix[,-1])
ARI_fuzzy.1[i]<-ARI.F(dat_surv$pop[dat_surv$J==1],prob_matrix[dat_surv$J==1,-1])
ARI_fuzzy.2[i]<-ARI.F(dat_surv$pop[dat_surv$J==2],prob_matrix[dat_surv$J==2,-1])
#hard
prob_df <- zeta_DP %>%
  mutate(item = row_number()) %>%
  pivot_longer(-item, values_to = "class") %>%
  count(item, class) %>%
  group_by(item) %>%
  mutate(prob = n / sum(n)) %>%
  select(-n)
hard_prob <- prob_df %>%
  group_by(item) %>%
  mutate(prob = if_else(prob == max(prob), 1, 0)) %>%
  ungroup()
hard_matrix <- hard_prob %>%
  pivot_wider(names_from = class, values_from = prob, values_fill = 0)
ARI_Hard.all[i]<-ARI.F(dat_surv$pop,hard_matrix[,-1])
ARI_Hard.1[i]<-ARI.F(dat_surv$pop[dat_surv$J==1],hard_matrix[dat_surv$J==1,-1])
ARI_Hard.2[i]<-ARI.F(dat_surv$pop[dat_surv$J==2],hard_matrix[dat_surv$J==2,-1])
}

ARI_DP<-data.frame(cbind(ARI_fuzzy.all,ARI_fuzzy.1,ARI_fuzzy.2,
                          ARI_Hard.all,ARI_Hard.1,ARI_Hard.2))

##############################DP_each
ARI_fuzzy.all<-rep(0,100)
ARI_fuzzy.1<-rep(0,100)
ARI_fuzzy.2<-rep(0,100)
ARI_Hard.all<-rep(0,100)
ARI_Hard.1<-rep(0,100)
ARI_Hard.2<-rep(0,100)
for (i in 1:100) {
load(paste0("DP_Each1_sim_Comp",i,".rdata"))
zeta_DPeach1<-as.data.frame(t(jags.m$BUGSoutput$sims.list$zeta[1:2000,]))
load(paste0("DP_Each2_sim_Comp",i,".rdata"))
zeta_DPeach2<-as.data.frame(t(jags.m$BUGSoutput$sims.list$zeta[1:2000,]))
zeta_DPeach<-rbind(zeta_DPeach1,zeta_DPeach2)
load(paste0("/Users/bohanjia/Desktop/research/gsr/Sim_Final/sim_100_train/dat_train",i,".rdata"))
dat_lmm<-dat%>%select(-event_time,-status)
dat_surv<-dat%>%group_by(id)%>%slice(1)%>%ungroup()%>%
    select(-y,-t)
df_long <- zeta_DPeach %>%
  mutate(item = row_number()) %>%
  pivot_longer(-item, names_to = "sample", values_to = "class")
prob_df <- df_long %>%
  group_by(item, class) %>%
  summarise(prob = n() / 2000, .groups = "drop")
prob_matrix <- prob_df %>%
  pivot_wider(names_from = class, values_from = prob, values_fill = 0)
#fuzzy
ARI_fuzzy.1[i]<-ARI.F(dat_surv$pop[dat_surv$J==1],prob_matrix[dat_surv$J==1,-1])
ARI_fuzzy.2[i]<-ARI.F(dat_surv$pop[dat_surv$J==2],prob_matrix[dat_surv$J==2,-1])
ARI_fuzzy.all[i]<-(ARI_fuzzy.1[i]*nrow(zeta_DPeach1)*(nrow(zeta_DPeach1)-1)/2+ 
                    ARI_fuzzy.2[i]*nrow(zeta_DPeach2)*(nrow(zeta_DPeach2)-1)/2)/(
                      nrow(zeta_DPeach1)*(nrow(zeta_DPeach1)-1)/2+
                        nrow(zeta_DPeach2)*(nrow(zeta_DPeach2)-1)/2 )
#hard
prob_df <- zeta_DPeach %>%
  mutate(item = row_number()) %>%
  pivot_longer(-item, values_to = "class") %>%
  count(item, class) %>%
  group_by(item) %>%
  mutate(prob = n / sum(n)) %>%
  select(-n)
hard_prob <- prob_df %>%
  group_by(item) %>%
  mutate(prob = if_else(prob == max(prob), 1, 0)) %>%
  ungroup()
hard_matrix <- hard_prob %>%
  pivot_wider(names_from = class, values_from = prob, values_fill = 0)
ARI_Hard.1[i]<-ARI.F(dat_surv$pop[dat_surv$J==1],hard_matrix[dat_surv$J==1,-1])
ARI_Hard.2[i]<-ARI.F(dat_surv$pop[dat_surv$J==2],hard_matrix[dat_surv$J==2,-1])
ARI_Hard.all[i]<-(ARI_Hard.1[i]*nrow(zeta_DPeach1)*(nrow(zeta_DPeach1)-1)/2+ 
  ARI_Hard.2[i]*nrow(zeta_DPeach2)*(nrow(zeta_DPeach2)-1)/2)/(
    nrow(zeta_DPeach1)*(nrow(zeta_DPeach1)-1)/2+
      nrow(zeta_DPeach2)*(nrow(zeta_DPeach2)-1)/2 )
}
ARI_DP_Each<-data.frame(cbind(ARI_fuzzy.all,ARI_fuzzy.1,ARI_fuzzy.2,
                              ARI_Hard.all,ARI_Hard.1,ARI_Hard.2))
colMeans(ARI_DP_Each)
colMeans(ARI_HDP)
colMeans(ARI_DP)

save(ARI_HDP, file = paste0("ARI_HDP.rdata"))
save(ARI_DP, file = paste0("ARI_DP.rdata"))
save(ARI_DP_Each, file = paste0("ARI_DP_Each.rdata"))

load("ARI_HDP.rdata")
load("ARI_DP.rdata")
load("ARI_DP_Each.rdata")
#ARI
colMeans(ARI_HDP)
colMeans(ARI_DP_Each)
colMeans(ARI_DP)
#whole pop
sd(ARI_HDP$ARI_Hard.all)
sd(ARI_DP_Each$ARI_Hard.all)
sd(ARI_DP$ARI_Hard.all)
wilcox.test(ARI_HDP$ARI_Hard.all, ARI_DP_Each$ARI_Hard.all,
            paired = TRUE, exact = FALSE)
wilcox.test(ARI_HDP$ARI_Hard.all, ARI_DP$ARI_Hard.all,
            paired = TRUE, exact = FALSE)


#whole pop
sd(ARI_HDP$ARI_fuzzy.all)
sd(ARI_DP_Each$ARI_fuzzy.all)
sd(ARI_DP$ARI_fuzzy.all)
wilcox.test(ARI_HDP$ARI_fuzzy.all, ARI_DP_Each$ARI_fuzzy.all,
            paired = TRUE, exact = FALSE)
wilcox.test(ARI_HDP$ARI_fuzzy.all, ARI_DP$ARI_fuzzy.all,
            paired = TRUE, exact = FALSE)
