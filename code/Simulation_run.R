# Simulated data

library("rgpdd")
library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")
library("purrr")

source("./code/GPDD_stability_functions.R")

sims=read.csv("./data/ChaosMetaAnalysisSimulatedData.csv")
sims_d=select(sims, ID, Model, TimeStep, NoiseLevel, Classification, Sim.1) %>% 
  group_by(ID,Model,Classification,NoiseLevel) %>%  mutate(TSlength=length(Sim.1)) %>% ungroup() %>% 
  group_by(ID,Model,Classification,NoiseLevel,TSlength) %>% nest() %>% 
  mutate(data=map(data, as.data.frame))
modelorder=unique(arrange(sims_d, Classification, Model)$Model)
sims_d$Model=factor(sims_d$Model, levels=modelorder)

#regression method
gpdd_results=select(sims_d, ID)
sims_results$regLE=map(sims_d$data, regLE, y="Sim.1")
sims_d$LEreg=map_dbl(sims_results$regLE, ~.x$LEreg)
sims_d$LEreg_se=map_dbl(sims_results$regLE, ~.x$LEreg_se)

ggplot(sims_d, aes(x=NoiseLevel, y=LEreg, color=Classification)) +
  facet_grid(TSlength~.) + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02)) + theme_bw()

#smap (model 1)
sims_results$modelresults1=map(sims_d$data, smap_model, y="Sim.1", ylog=F)
sims_results$jacobians1=map(sims_results$modelresults1, getJacobians)
sims_results$stability1=map(sims_results$jacobians1, getStability)

sims_d$gle=map_dbl(sims_results$stability1, ~.x$gle)
sims_d$Ebest=map_dbl(sims_results$modelresults1, ~.x$modelstats$E)
sims_d$thetabest=map_dbl(sims_results$modelresults1, ~.x$modelstats$theta)
sims_d$R1a=map_dbl(sims_results$modelresults1, ~.x$modelstats$R2abund)

ggplot(sims_d, aes(x=NoiseLevel, y=gle, color=Classification)) +
  facet_grid(TSlength~., scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02)) + theme_bw()

ggplot(sims_d, aes(x=Model, y=gle, color=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02)) + theme_bw() + xlabvert

