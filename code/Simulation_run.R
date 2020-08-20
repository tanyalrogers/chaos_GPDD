# Simulated data

library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")
library("purrr")

source("./code/GPDD_stability_functions.R")
#load(file = "./data/sims_results.Rdata")
source("~/GRAD SCHOOL/R reference/ggplot themes rogers.R")


sims=read.csv("./data/ChaosMetaAnalysisSimulatedDataCORRECTED.csv")
sims$Classification=recode(sims$Classification, Periodic="periodic")
sims_d=select(sims, ID, Model, TimeStep, NoiseLevel, Classification, Sim.1:Sim.2) %>% 
  gather(SimNumber, Value, Sim.1:Sim.2) %>% 
  group_by(ID,Model,Classification,NoiseLevel,SimNumber) %>%  mutate(TSlength=length(Value)) %>% ungroup() %>% 
  group_by(ID,Model,Classification,NoiseLevel,TSlength,SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame))
modelorder=unique(arrange(sims_d, Classification, Model)$Model)
sims_d$Model=factor(sims_d$Model, levels=modelorder)

#plot timeseries
sims_plot=filter(sims_d, SimNumber=="Sim.1" & TSlength==50 & NoiseLevel==0.1)
par(mfrow=c(2,3))
for(i in 1:nrow(sims_plot)) {
  dtemp=sims_plot$data[[i]]
  plot(Value~TimeStep, data=dtemp, main=sims_plot$Model[i], type="l")
  #pacf(dtemp$Value,main=sims_plot$Model[i])
}

sims_results=select(sims_d, ID, SimNumber)

#best hyperparameters
sims_results$hyperpars=map(sims_d$data, besthyper, y="Value")
sims_d=cbind(sims_d, bind_rows(sims_results$hyperpars))

ggplot(sims_d, aes(x=Model, y=bestE, color=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point() + #geom_boxplot() + 
  theme_bw() + xlabvert
ggplot(sims_d, aes(x=Model, y=bestTau, color=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point() + #geom_boxplot() + 
  theme_bw() + xlabvert

#regression method
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

#confidence interval on mean lle
sims_d$lle_avg=map_dbl(sims_results$stability1, ~.x$lle_avg)
sims_d$lle_ci90lower=map_dbl(sims_results$stability1, function(x) {
  std=sd(x$lle, na.rm=T)
  n=length(which(!is.na(x$lle)))
  ci=x$lle_avg+std/sqrt(n)*qt(p=0.05, df=n-1)})
sims_d$lle_ci90upper=map_dbl(sims_results$stability1, function(x) {
  std=sd(x$lle, na.rm=T)
  n=length(which(!is.na(x$lle)))
  ci=x$lle_avg+std/sqrt(n)*qt(p=0.95, df=n-1)})
sims_d$lle_class=ifelse(sims_d$lle_ci90lower>0, "chaotic", "not chaotic")


#
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
  group_by(Classification,NoiseLevel,TSlength) %>% 
  summarize(gle_pp=length(which(gle>0))/length(gle),
            gle_pp.01=length(which(gle>0.01))/length(gle),
            gle_pp.05=length(which(gle>0.05))/length(gle),
            reg_pp=length(which(LEreg>0))/length(LEreg),
            lle_avg_pp=length(which(lle_avg>0))/length(lle_avg),
            lle_pp=length(which(lle_class=="chaotic"))/length(lle_class))

ggplot(sims_summary, aes(x=Model, y=lle_pp, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + 
  geom_bar(stat = "identity") + theme_bw() + xlabvert
ggplot(sims_summary, aes(x=Model, y=lle_avg_pp, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + 
  geom_bar(stat = "identity") + theme_bw() + xlabvert
ggplot(sims_summary, aes(x=Model, y=gle_pp, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + 
  geom_bar(stat = "identity") + theme_bw() + xlabvert
ggplot(sims_summary, aes(x=Model, y=gle_pp.05, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + 
  geom_bar(stat = "identity") + theme_bw() + xlabvert
ggplot(sims_summary, aes(x=Model, y=reg_pp, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_bar(position = position_dodge2(), stat = "identity") + theme_bw() + xlabvert

ggplot(sims_summary, aes(x=NoiseLevel, y=gle_pp, fill=Classification)) +
  facet_grid(TSlength~Classification) + 
  geom_bar(stat = "identity", position = position_dodge2(), show.legend = F) + theme_bw() + xlabvert
ggplot(sims_summary, aes(x=NoiseLevel, y=lle_pp, fill=Classification)) +
  facet_grid(TSlength~Classification) + 
  geom_bar(stat = "identity", position = position_dodge2(), show.legend = F) + theme_bw() + xlabvert
ggplot(sims_summary, aes(x=NoiseLevel, y=lle_avg_pp, fill=Classification)) +
  facet_grid(TSlength~Classification) + 
  geom_bar(stat = "identity", position = position_dodge2(), show.legend = F) + theme_bw() + xlabvert
