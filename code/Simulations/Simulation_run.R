# Applies Jacobian and Direct LE chaos detection methods to simulated data
# Tanya Rogers

library(dplyr)
library(tidyr)
library(rEDM)
library(purrr)
library(furrr)
library(ggplot2)

source("./code/Methods/LE_ChaosDetectionMethods.R")

#to run in parallel
plan(multisession, workers = 4)

#### Simulations with known E and tau (embedding dataset) ####

sims=read.csv("./data/simulation_dataset_embedding.csv")
sims_d=select(sims, TSlength=Time.Series.Length, Known.Tau, Known.E, TimeStep=Time.Step, Sim.1:Sim.10) %>% 
  gather(SimNumber, Value, Sim.1:Sim.10) %>% 
  group_by(TSlength,Known.Tau, Known.E,SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame)) %>% arrange(Known.E, Known.Tau)
sims_d$modelresults1=map(sims_d$data, smap_model_options, y="Value", model=1)
sims_d$modelresults2=map(sims_d$data, smap_model_options, y="Value", model=2)
sims_d$modelresults4=map(sims_d$data, smap_model_options, y="Value", model=4)
sims_d$R1a=map_dbl(sims_d$modelresults1, ~.x$modelstats$R2abund)
sims_d$Ebest1=map_dbl(sims_d$modelresults1, ~.x$modelstats$E)
sims_d$taubest1=map_dbl(sims_d$modelresults1, ~.x$modelstats$tau)
sims_d$R2a=map_dbl(sims_d$modelresults2, ~.x$modelstats$R2abund)
sims_d$Ebest2=map_dbl(sims_d$modelresults2, ~.x$modelstats$E)
sims_d$taubest2=map_dbl(sims_d$modelresults2, ~.x$modelstats$tau)
sims_d$R4a=map_dbl(sims_d$modelresults4, ~.x$modelstats$R2abund)
sims_d$Ebest4=map_dbl(sims_d$modelresults4, ~.x$modelstats$E)
sims_d$taubest4=map_dbl(sims_d$modelresults4, ~.x$modelstats$tau)
sims_d$thetabest4=map_dbl(sims_d$modelresults4, ~.x$modelstats$theta)
sims_d$bestmodel=select(sims_d,R1a,R2a,R4a) %>% apply(1,which.max)

ggplot(sims_d, aes(x=Ebest4, y=taubest4)) +
  facet_grid(Known.E~Known.Tau) + 
  geom_point(aes(x=Known.E, y=Known.Tau), color="red", size=3) +
  geom_point(alpha=0.5) + theme_bw()

# save(sims_d,  file = "./data/sims_results_knownEtau.Rdata")

#### Test dataset with known dynamics ####

sims=read.csv("./data/simulation_dataset_test.csv")

sims_d=select(sims, ID, Model, TimeStep, NoiseLevel, Classification, Sim.1:Sim.100) %>% 
  gather(SimNumber, Value, Sim.1:Sim.100) %>% 
  group_by(ID,Model,Classification,NoiseLevel,SimNumber) %>%  mutate(TSlength=length(Value)) %>% ungroup() %>% 
  group_by(ID,Model,Classification,NoiseLevel,TSlength,SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame)) %>% ungroup()
modelorder=unique(arrange(sims_d, Classification, Model)$Model)
sims_d$Model=factor(sims_d$Model, levels=modelorder)

#set up results df
sims_results=select(sims_d, ID, SimNumber, Model)

# Jacobian LE method

#get hyperparameters (this takes a long time)
sims_results$hpar1=future_map(sims_d$data, besthyper, y="Value", ylog=F, pgr="none")
#get model output
sims_results$modelresults1=map2(sims_d$data, sims_results$hpar1, smap_model_options, y="Value", model=1)
#get R2
sims_d$R1a=map_dbl(sims_results$modelresults1, ~.x$modelstats$R2abund)

#fit alternate models for strictly positive data
# #1 and 3 are same, 2 and 4 are same
logmodels=c("PredatorPreyPeriodic", "PredatorPreyChaotic", "HostParParPeriodic", "HostParParChaotic")
sims_log=filter(sims_d, Model %in% logmodels)
sims_log_results=filter(sims_results, Model %in% logmodels)
sims_log_results$hpar4=future_map(sims_log$data, besthyper, y="Value", ylog=F, pgr="gr")
sims_log_results$hpar5=future_map(sims_log$data, besthyper, y="Value", ylog=T, pgr="gr")
sims_log_results$modelresults4=map2(sims_log$data, sims_log_results$hpar4, smap_model_options, y="Value", model=4)
sims_log_results$modelresults5=map2(sims_log$data, sims_log_results$hpar5, smap_model_options, y="Value", model=5)
sims_log$R1a=map_dbl(sims_log_results$modelresults1, ~.x$modelstats$R2abund)
sims_log$R4a=map_dbl(sims_log_results$modelresults4, ~.x$modelstats$R2abund)
sims_log$R5a=map_dbl(sims_log_results$modelresults5, ~.x$modelstats$R2abund)
sims_log$bestR2=select(sims_log,R1a,R4a,R5a) %>% apply(1,max)
sims_log$bestmodel=select(sims_log,R1a,R4a,R5a) %>% apply(1,which.max)
sims_log_results$bestmodel=select(sims_log,R1a,R4a,R5a) %>% apply(1,which.max)
sims_log_results$modelresultsbest=cbind(select(sims_log_results, modelresults3, modelresults4, modelresults5),sims_log$bestmodel) %>% apply(1, function(x) {m=as.numeric(x["sims_log$bestmodel"]); x[m][[1]]})

#join log models with others
sims_results=left_join(sims_results, select(sims_log_results, ID, SimNumber, modelresultsbest, bestmodel), by=c("ID", "SimNumber"))
sims_results$modelresultsbest=ifelse(sims_d$Model %in% logmodels,
                                     sims_results$modelresultsbest, 
                                     sims_results$modelresults1)

#get stability
sims_results$jacobians=map(sims_results$modelresultsbest, getJacobians)
sims_results$LEshift=map2(sims_results$modelresultsbest, sims_results$jacobians, LEshift)

#pull results
sims_d$Ebest=map_dbl(sims_results$modelresultsbest, ~.x$modelstats$E)
sims_d$taubest=map_dbl(sims_results$modelresultsbest, ~.x$modelstats$tau)
sims_d$thetabest=map_dbl(sims_results$modelresultsbest, ~.x$modelstats$theta)
sims_d$R2best=map_dbl(sims_results$modelresultsbest, ~.x$modelstats$R2abund)
sims_d$gle=map_dbl(sims_results$stability, ~.x$gle)
sims_d$minmean=map_dbl(sims_results$LEshift, ~.x$minmean)
sims_d$minci=map_dbl(sims_results$LEshift, ~.x$minci)

# Direct LE method

sims_results$regLE=map2(sims_d$data, sims_results$modelresultsbest, regLE, y="Value")
sims_d$LEreg=map_dbl(sims_results$regLE, ~.x$LEreg)
sims_d$LEreg_se=map_dbl(sims_results$regLE, ~.x$LEreg_se)

#save(sims_d, sims_results, file = "./data/sims_results_update2.Rdata")
#save(sims_log, sims_log_results, sims2_log, sims2_log_results, file = "./data/sims_results_log_update.Rdata")

# Export results

#get best model form
sims_d$modelform=map_chr(sims_results$modelresultsbest, ~.x$form)
#reclass noise level for stochastic ts
sims_d$NoiseLevel2=ifelse(sims_d$NoiseLevel==0, 0.01, sims_d$NoiseLevel)
sims_d$Classification2=ifelse(sims_d$Classification=="chaotic", "chaotic", "not chaotic")
#class LEs Direct (regression) method
sims_d$LEregmin=sims_d$LEreg-1.96*sims_d$LEreg_se
sims_d$LEregclass=ifelse(sims_d$LEregmin>0.01, "chaotic", "not chaotic")
#class LEs Jacobian method
sims_d$LEclass=ifelse(sims_d$minci>0.01, "chaotic", "not chaotic")

#export E and tau for other analyses
Eexport=spread(select(sims_d, ID:SimNumber, Ebest), SimNumber, Ebest) %>% 
  select(ID:TSlength, paste0("Sim.",1:100))
tauexport=spread(select(sims_d, ID:SimNumber, taubest), SimNumber, taubest) %>% 
  select(ID:TSlength, paste0("Sim.",1:100))
write.csv(Eexport,"./data/sims_test_E.csv", row.names = F)
write.csv(tauexport,"./data/sims_test_tau.csv", row.names = F)

#export results
dexport=select(sims_d, ID:SimNumber, E=Ebest, tau=taubest, theta=thetabest, R2=R2best, modelform, 
               LEmean=minmean, LEmin=minci, LEreg, LEreg_se, NoiseLevel2:LEclass)
write.csv(dexport,"./data/sims_test_results.csv", row.names = F)

#### Validation dataset with known dynamics ####

sims=read.csv("./data/simulation_dataset_validation.csv")

sims_v=select(sims, ID, Model, TimeStep, NoiseLevel, TSlength=TimeSeriesLength, Classification, Sim.1:Sim.100) %>% 
  gather(SimNumber, Value, Sim.1:Sim.100) %>% 
  group_by(ID,Model,Classification,NoiseLevel,TSlength,SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame)) %>% ungroup()
modelorder=unique(arrange(sims_v, Classification, Model)$Model)
sims_v$Model=factor(sims_v$Model, levels=modelorder)

# Jacobian LE method

#set up results df
sims_vresults=select(sims_v, ID, SimNumber, Model)
#get hyperparameters (this takes a long time)
sims_vresults$hpar1=future_map(sims_v$data, besthyper, y="Value", ylog=F, pgr="none")
#get model output
sims_vresults$modelresults1=map2(sims_v$data, sims_vresults$hpar1, smap_model_options, y="Value", model=1)
#get R2
sims_v$R1a=map_dbl(sims_vresults$modelresults1, ~.x$modelstats$R2abund)

#fit alt models for strictly positive data
# #1 and 3 are same, 2 and 4 are same
logmodels=c("CompetitionChaotic", "CompetitionPeriodic", "TinkerbellChaotic", "TinkerbellPeriodic")
sims_vlog=filter(sims_v, Model %in% logmodels)
sims_vlog_results=filter(sims_vresults, Model %in% logmodels)
sims_vlog_results$hpar4=future_map(sims_vlog$data, besthyper, y="Value", ylog=F, pgr="gr")
sims_vlog_results$hpar5=future_map(sims_vlog$data, besthyper, y="Value", ylog=T, pgr="gr")
sims_vlog_results$modelresults4=map2(sims_vlog$data, sims_vlog_results$hpar4, smap_model_options, y="Value", model=4)
sims_vlog_results$modelresults5=map2(sims_vlog$data, sims_vlog_results$hpar5, smap_model_options, y="Value", model=5)
sims_vlog$R1a=map_dbl(sims_vlog_results$modelresults1, ~.x$modelstats$R2abund)
sims_vlog$R4a=map_dbl(sims_vlog_results$modelresults4, ~.x$modelstats$R2abund)
sims_vlog$R5a=map_dbl(sims_vlog_results$modelresults5, ~.x$modelstats$R2abund)
sims_vlog$bestR2=select(sims_vlog,R1a,R4a,R5a) %>% apply(1,max)
sims_vlog$bestmodel=select(sims_vlog,R1a,R4a,R5a) %>% apply(1,which.max)
sims_vlog_results$bestmodel=select(sims_vlog,R1a,R4a,R5a) %>% apply(1,which.max)
sims_vlog_results$modelresultsbest=cbind(select(sims_vlog_results, modelresults1, modelresults4, modelresults5),sims_vlog$bestmodel) %>% apply(1, function(x) {m=as.numeric(x["sims_vlog$bestmodel"]); x[m][[1]]})

#join log models with others
sims_vresults=left_join(sims_vresults, select(sims_vlog_results, ID, SimNumber, modelresultsbest, bestmodel), by=c("ID", "SimNumber"))
sims_vresults$modelresultsbest=ifelse(sims_v$Model %in% logmodels,
                                      sims_vresults$modelresultsbest, 
                                      sims_vresults$modelresults1)

#get stability
sims_vresults$jacobians=map(sims_vresults$modelresultsbest, getJacobians)
sims_vresults$LEshift=map2(sims_vresults$modelresultsbest, sims_vresults$jacobians, LEshift)

#pull results
sims_v$Ebest=map_dbl(sims_vresults$modelresultsbest, ~.x$modelstats$E)
sims_v$taubest=map_dbl(sims_vresults$modelresultsbest, ~.x$modelstats$tau)
sims_v$thetabest=map_dbl(sims_vresults$modelresultsbest, ~.x$modelstats$theta)
sims_v$R2best=map_dbl(sims_vresults$modelresultsbest, ~.x$modelstats$R2abund)
sims_v$minmean=map_dbl(sims_vresults$LEshift, ~.x$minmean)
sims_v$minci=map_dbl(sims_vresults$LEshift, ~.x$minci)

# Direct LE method

sims_vresults$regLE=map2(sims_v$data, sims_vresults$modelresultsbest, regLE, y="Value")
sims_v$LEreg=map_dbl(sims_vresults$regLE, ~.x$LEreg)
sims_v$LEreg_se=map_dbl(sims_vresults$regLE, ~.x$LEreg_se)

#save(sims_v, sims_vresults, file = "./data/sims_validation2.Rdata")
#save(sims_vlog, sims_vlog_results, file = "./data/sims_validation_log.Rdata")

# Export results

#get best model form
sims_v$modelform=map_chr(sims_vresults$modelresultsbest, ~.x$form)
#reclass noise level for stochastic ts
sims_v$NoiseLevel2=ifelse(sims_v$NoiseLevel==0, 0.01, sims_v$NoiseLevel)
sims_v$Classification2=ifelse(sims_v$Classification=="chaotic", "chaotic", "not chaotic")
#class LEs regression method
sims_v$LEregmin=sims_v$LEreg-1.96*sims_v$LEreg_se
sims_v$LEregclass=ifelse(sims_v$LEregmin>0.01, "chaotic", "not chaotic")
#class LEs jacobian method
sims_v$LEclass=ifelse(sims_v$minci>0.01, "chaotic", "not chaotic")

#export E and tau for other analyses
Eexport=spread(select(sims_v, ID:SimNumber, Ebest), SimNumber, Ebest) %>% 
  select(ID:TSlength, paste0("Sim.",1:100))
tauexport=spread(select(sims_v, ID:SimNumber, taubest), SimNumber, taubest) %>% 
  select(ID:TSlength, paste0("Sim.",1:100))
write.csv(Eexport,"./data/sims_validation_E.csv", row.names = F)
write.csv(tauexport,"./data/sims_validation_tau.csv", row.names = F)

vexport=select(sims_v, ID:SimNumber, E=Ebest, tau=taubest, theta=thetabest, R2=R2best, modelform, 
               LEmean=minmean, LEmin=minci, LEreg, LEreg_se, NoiseLevel2:LEclass)
write.csv(vexport,"./data/sims_validation_results.csv", row.names = F)

#### Observation and process noise dataset with known dynamics ####

sims=read.csv("./data/simulation_dataset_noise_test.csv")

sims_n=select(sims, ID, Model, TimeStep, ObsNoise, ProcessNoise, TSlength=TimeSeriesLength, Classification, Sim.1:Sim.100) %>% 
  gather(SimNumber, Value, Sim.1:Sim.100) %>% 
  group_by(ID,Model,Classification,ObsNoise, ProcessNoise,TSlength,SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame)) %>% ungroup()
modelorder=unique(arrange(sims_n, Classification, Model)$Model)
sims_n$Model=factor(sims_n$Model, levels=modelorder)

# Jacobian LE method

#set up results df
sims_nresults=select(sims_n, ID, SimNumber, Model)
#get hyperparameters (this takes a long time)
starttime=Sys.time()
sims_nresults$hpar1=future_map(sims_n$data, besthyper, y="Value", ylog=F, pgr="none")
endtime=Sys.time(); endtime-starttime
starttime=Sys.time()
sims_nresults$hpar4=future_map(sims_n$data, besthyper, y="Value", ylog=F, pgr="gr")
sims_nresults$hpar5=future_map(sims_n$data, besthyper, y="Value", ylog=T, pgr="gr")
endtime=Sys.time(); endtime-starttime
#get model output
sims_nresults$modelresults1=map2(sims_n$data, sims_nresults$hpar1, smap_model_options, y="Value", model=1)
sims_nresults$modelresults4=map2(sims_n$data, sims_nresults$hpar4, smap_model_options, y="Value", model=4)
sims_nresults$modelresults5=map2(sims_n$data, sims_nresults$hpar5, smap_model_options, y="Value", model=5)
#get R2
sims_n$R1a=map_dbl(sims_nresults$modelresults1, ~.x$modelstats$R2abund)
sims_n$R4a=map_dbl(sims_nresults$modelresults4, ~.x$modelstats$R2abund)
sims_n$R5a=map_dbl(sims_nresults$modelresults5, ~.x$modelstats$R2abund)

sims_n$bestR2=select(sims_n,R1a,R4a,R5a) %>% apply(1,max)
sims_n$bestmodel=select(sims_n,R1a,R4a,R5a) %>% apply(1,which.max)
sims_nresults$bestmodel=select(sims_n,R1a,R4a,R5a) %>% apply(1,which.max)
sims_nresults$modelresultsbest=cbind(select(sims_nresults, modelresults1, modelresults4, modelresults5),sims_n$bestmodel) %>% apply(1, function(x) {m=as.numeric(x["sims_n$bestmodel"]); x[m][[1]]})

#get stability
sims_nresults$jacobians=map(sims_nresults$modelresultsbest, getJacobians)
sims_nresults$LEshift=map2(sims_nresults$modelresultsbest, sims_nresults$jacobians, LEshift)

#pull results
sims_n$Ebest=map_dbl(sims_nresults$modelresultsbest, ~.x$modelstats$E)
sims_n$taubest=map_dbl(sims_nresults$modelresultsbest, ~.x$modelstats$tau)
sims_n$thetabest=map_dbl(sims_nresults$modelresultsbest, ~.x$modelstats$theta)
sims_n$R2best=map_dbl(sims_nresults$modelresultsbest, ~.x$modelstats$R2abund)
sims_n$minmean=map_dbl(sims_nresults$LEshift, ~.x$minmean)
sims_n$minci=map_dbl(sims_nresults$LEshift, ~.x$minci)

# Direct LE method

sims_nresults$regLE=map2(sims_n$data, sims_nresults$modelresultsbest, regLE, y="Value")
sims_n$LEreg=map_dbl(sims_nresults$regLE, ~.x$LEreg)
sims_n$LEreg_se=map_dbl(sims_nresults$regLE, ~.x$LEreg_se)

save(sims_n, sims_nresults, file = "./data/sims_noise.Rdata")

# Export results

#get best model form
sims_n$modelform=map_chr(sims_nresults$modelresultsbest, ~.x$form)
sims_n$Classification2=ifelse(sims_n$Classification=="chaotic", "chaotic", "not chaotic")
#class LEs regression method
sims_n$LEregmin=sims_n$LEreg-1.96*sims_n$LEreg_se
sims_n$LEregclass=ifelse(sims_n$LEregmin>0.01, "chaotic", "not chaotic")
#class LEs jacobian method
sims_n$LEclass=ifelse(sims_n$minci>0.01, "chaotic", "not chaotic")

#export E and tau for other analyses
Eexport=spread(select(sims_n, ID:SimNumber, Ebest), SimNumber, Ebest) %>% 
  select(ID:TSlength, paste0("Sim.",1:100))
tauexport=spread(select(sims_n, ID:SimNumber, taubest), SimNumber, taubest) %>% 
  select(ID:TSlength, paste0("Sim.",1:100))
write.csv(Eexport,"./data/sims_noise_E.csv", row.names = F)
write.csv(tauexport,"./data/sims_noise_tau.csv", row.names = F)

nexport=select(sims_n, ID:SimNumber, E=Ebest, tau=taubest, theta=thetabest, R2=R2best, modelform, 
               LEmean=minmean, LEmin=minci, LEreg, LEreg_se, Classification2, LEregmin:LEclass)
write.csv(nexport,"./data/sims_noise_results.csv", row.names = F)
