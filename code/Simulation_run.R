# Simulated data

library("rgpdd")
library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")
library("purrr")

source("./code/GPDD_stability_functions.R")

sims=read.csv("./data/ChaosMetaAnalysisSimulatedData.csv")
sims_d=select(sims, ID, Model, TimeStep, NoiseLevel, Classification, Sim.1:Sim.10) %>% 
  gather(SimNumber, Value, Sim.1:Sim.10) %>% 
  group_by(ID,Model,Classification,NoiseLevel,SimNumber) %>%  mutate(TSlength=length(Value)) %>% ungroup() %>% 
  group_by(ID,Model,Classification,NoiseLevel,TSlength,SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame))
modelorder=unique(arrange(sims_d, Classification, Model)$Model)
sims_d$Model=factor(sims_d$Model, levels=modelorder)

#regression method
sims_results=select(sims_d, ID, SimNumber)
sims_results$regLE=map(sims_d$data, regLE, y="Value")
sims_d$LEreg=map_dbl(sims_results$regLE, ~.x$LEreg)
sims_d$LEreg_se=map_dbl(sims_results$regLE, ~.x$LEreg_se)

ggplot(sims_d, aes(x=NoiseLevel, y=LEreg, color=Classification)) +
  facet_grid(TSlength~.) + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02)) + theme_bw()

ggplot(sims_d, aes(x=Model, y=LEreg, color=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point() + #geom_boxplot() + 
  theme_bw() + xlabvert

#smap (model 1)
sims_results$modelresults1=map(sims_d$data, smap_model, y="Value", ylog=F)
sims_results$jacobians1=map(sims_results$modelresults1, getJacobians)
sims_results$stability1=map(sims_results$jacobians1, getStability)

sims_d$gle=map_dbl(sims_results$stability1, ~.x$gle)
sims_d$Ebest=map_dbl(sims_results$modelresults1, ~.x$modelstats$E)
sims_d$thetabest=map_dbl(sims_results$modelresults1, ~.x$modelstats$theta)
sims_d$R1a=map_dbl(sims_results$modelresults1, ~.x$modelstats$R2abund)

save(sims_d, sims_results, file = "./data/sims_results.Rdata")

ggplot(sims_d, aes(x=NoiseLevel, y=gle, color=Classification)) +
  facet_grid(TSlength~., scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.05), alpha=0.4) + 
  theme_bw() + xlab("Noise Level") + legalpha
ggplot(sims_d, aes(x=factor(NoiseLevel), y=gle, fill=Classification)) +
  facet_grid(TSlength~., scales = "free_y") + geom_hline(yintercept = 0) +
  geom_violin(position = position_dodge(0.9), scale = "width") + 
  theme_bw() + xlab("Noise Level")

ggplot(sims_d, aes(x=Model, y=gle, color=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02)) + theme_bw() + xlabvert

#R squared, E, theta
ggplot(sims_d, aes(x=Model, y=R1a, color=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02), alpha=0.3) + theme_bw() + xlabvert + legalpha + ylab("R-squared")
ggplot(sims_d, aes(x=Model, y=Ebest, color=Classification)) +
  facet_grid(TSlength~NoiseLevel) + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02), alpha=0.3) + theme_bw() + xlabvert + legalpha
ggplot(sims_d, aes(x=Model, y=thetabest, color=Classification)) +
  facet_grid(TSlength~NoiseLevel) + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02), alpha=0.3) + theme_bw() + xlabvert + legalpha

#histograms
tslengths=unique(sims_d$TSlength)
for(i in 1:length(tslengths)) {
  print(ggplot(filter(sims_d, TSlength==tslengths[i]), aes(x=gle, fill=Classification)) +
    facet_grid(Classification~NoiseLevel, scales = "free_y") + geom_hline(yintercept = 0) +
    geom_histogram(boundary = 0, binwidth = 0.1, show.legend = F) + geom_vline(xintercept = 0) +
    theme_bw() + ggtitle(paste("TSlength =", tslengths[i])) + xlim(c(-1,1)))
}
ggplot(sims_d, aes(x=gle, color=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_freqpoly(boundary = 0, binwidth = 0.2, size=1) + geom_vline(xintercept = 0) +
  theme_bw() + xlim(c(-1,1))


sims_summary=sims_d %>% select(-data) %>% 
  group_by(ID,Model,Classification,NoiseLevel,TSlength) %>% 
  summarize(gle_pp=length(which(gle>0))/length(gle),
            gle_pp.01=length(which(gle>0.01))/length(gle),
            gle_pp.05=length(which(gle>0.05))/length(gle),
            reg_pp=length(which(LEreg>0))/length(LEreg))

ggplot(sims_summary, aes(x=Model, y=gle_pp, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + 
  geom_bar(stat = "identity") + theme_bw() + xlabvert
ggplot(sims_summary, aes(x=Model, y=gle_pp.05, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + 
  geom_bar(stat = "identity") + theme_bw() + xlabvert

ggplot(sims_summary, aes(x=Model, y=reg_pp, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_bar(position = position_dodge2(), stat = "identity") + theme_bw() + xlabvert
