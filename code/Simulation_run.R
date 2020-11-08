# Simulated data

library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")
library("purrr")
library("furrr")

source("./code/GPDD_stability_functions.R")
source("./code/ggplot themes rogers.R")

#to run in parallel
plan(multisession, workers = 4)

#plot timeseries ####
sims_plot=filter(sims_d, SimNumber=="Sim.3" & TSlength==50 & NoiseLevel==0.3)
par(mfrow=c(2,3))
for(i in 1:nrow(sims_plot)) {
  dtemp=sims_plot$data[[i]]
  plot(Value~TimeStep, data=dtemp, main=sims_plot$Model[i], type="l")
  #pacf(dtemp$Value,main=sims_plot$Model[i])
}

#test with known E and tau ####
sims=read.csv("./data/KnownEandTau2.csv")
sims_d=select(sims, TSlength=Time.Series.Length, Known.Tau, Known.E, TimeStep=Time.Step, Sim.1:Sim.10) %>% 
  gather(SimNumber, Value, Sim.1:Sim.10) %>% 
  group_by(TSlength,Known.Tau, Known.E,SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame)) %>% arrange(Known.E, Known.Tau)
sims_d$modelresults1=map(sims_d$data, smap_model_options, y="Value", model=1)
sims_d$modelresults2=map(sims_d$data, smap_model_options, y="Value", model=2)
sims_d$modelresults4=map(sims_d$data, smap_model_options, y="Value", model=4)
#sims_d$modelresults4=pmap(list(sims_d$data, Efix=sims_d$Known.E, taufix=sims_d$Known.Tau), smap_model_options, y="Value", model=4)
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

save(sims_d,  file = "./data/sims_results_knownEtau.Rdata")

#test with known dynamics ####
sims=read.csv("./data/ChaosMetaAnalysisSimulatedDataCORRECTED3.csv")
sims$Classification=recode(sims$Classification, Periodic="periodic")

#first 20 reps
sims_d=select(sims, ID, Model, TimeStep, NoiseLevel, Classification, Sim.1:Sim.20) %>% 
  gather(SimNumber, Value, Sim.1:Sim.20) %>% 
  group_by(ID,Model,Classification,NoiseLevel,SimNumber) %>%  mutate(TSlength=length(Value)) %>% ungroup() %>% 
  group_by(ID,Model,Classification,NoiseLevel,TSlength,SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame)) %>% ungroup()
  #left_join(sims_Etau, by=c("ID", "SimNumber"))
modelorder=unique(arrange(sims_d, Classification, Model)$Model)
sims_d$Model=factor(sims_d$Model, levels=modelorder)

#set up results df
sims_results=select(sims_d, ID, SimNumber, Model)

#get hyperparameters (this takes a long time)
sims_results$hpar1=future_map(sims_d$data, besthyper, y="Value", ylog=F, pgr="none")
#get model output
sims_results$modelresults1=map2(sims_d$data, sims_results$hpar1, smap_model_options, y="Value", model=1)
#get R2
sims_d$R1a=map_dbl(sims_results$modelresults1, ~.x$modelstats$R2abund)

#fit alt models for strictly positive data
# #1 and 3 are same, 2 and 4 are same
logmodels=c("PredatorPreyPeriodic", "PredatorPreyChaotic", "HostParParPeriodic", "HostParParChaotic")
sims_log=filter(sims_d, Model %in% logmodels)
sims_log_results=filter(sims_results, Model %in% logmodels)
sims_log_results$hpar2=future_map(sims_log$data, besthyper, y="Value", ylog=T, pgr="none")
sims_log_results$hpar3=future_map(sims_log$data, besthyper, y="Value", ylog=F, pgr="fd")
sims_log_results$hpar4=future_map(sims_log$data, besthyper, y="Value", ylog=F, pgr="gr")
sims_log_results$hpar5=future_map(sims_log$data, besthyper, y="Value", ylog=T, pgr="gr")
sims_log_results$modelresults2=map2(sims_log$data, sims_log_results$hpar2, smap_model_options, y="Value", model=2)
sims_log_results$modelresults3=map2(sims_log$data, sims_log_results$hpar3, smap_model_options, y="Value", model=3)
sims_log_results$modelresults4=map2(sims_log$data, sims_log_results$hpar4, smap_model_options, y="Value", model=4)
sims_log_results$modelresults5=map2(sims_log$data, sims_log_results$hpar5, smap_model_options, y="Value", model=5)
sims_log$R1a=map_dbl(sims_log_results$modelresults1, ~.x$modelstats$R2abund)
sims_log$R2a=map_dbl(sims_log_results$modelresults2, ~.x$modelstats$R2abund)
sims_log$R3a=map_dbl(sims_log_results$modelresults3, ~.x$modelstats$R2abund)
sims_log$R4a=map_dbl(sims_log_results$modelresults4, ~.x$modelstats$R2abund)
sims_log$R5a=map_dbl(sims_log_results$modelresults5, ~.x$modelstats$R2abund)
sims_log$bestR2=select(sims_log,R3a,R4a,R5a) %>% apply(1,max)
sims_log$bestmodel=select(sims_log,R3a,R4a,R5a) %>% apply(1,which.max)
sims_log_results$bestmodel=select(sims_log,R3a,R4a,R5a) %>% apply(1,which.max)
sims_log_results$modelresultsbest=cbind(select(sims_log_results, modelresults3, modelresults4, modelresults5),sims_log$bestmodel) %>% apply(1, function(x) {m=as.numeric(x["sims_log$bestmodel"]); x[m][[1]]})

#join log models with others
sims_results=left_join(sims_results, select(sims_log_results, ID, SimNumber, modelresultsbest, bestmodel), by=c("ID", "SimNumber"))
sims_results$modelresultsbest=ifelse(sims_d$Model %in% logmodels,
                                     sims_results$modelresultsbest, 
                                     sims_results$modelresults1)

#get stability
sims_results$jacobians=map(sims_results$modelresultsbest, getJacobians)
sims_results$stability=map2(sims_results$modelresultsbest, sims_results$jacobians, getStability)
sims_results$LEshift=map2(sims_results$modelresultsbest, sims_results$jacobians, LEshift)

#pull results
sims_d$Ebest=map_dbl(sims_results$modelresultsbest, ~.x$modelstats$E)
sims_d$taubest=map_dbl(sims_results$modelresultsbest, ~.x$modelstats$tau)
sims_d$thetabest=map_dbl(sims_results$modelresultsbest, ~.x$modelstats$theta)
sims_d$R2best=map_dbl(sims_results$modelresultsbest, ~.x$modelstats$R2abund)
sims_d$gle=map_dbl(sims_results$stability, ~.x$gle)
sims_d$minmean=map_dbl(sims_results$LEshift, ~.x$minmean)
sims_d$minci=map_dbl(sims_results$LEshift, ~.x$minci)

#remaining reps ####
#21-50
sims2_d=select(sims, ID, Model, TimeStep, NoiseLevel, Classification, Sim.21:Sim.50) %>% 
  gather(SimNumber, Value, Sim.21:Sim.50) %>% 
  group_by(ID,Model,Classification,NoiseLevel,SimNumber) %>%  mutate(TSlength=length(Value)) %>% ungroup() %>% 
  group_by(ID,Model,Classification,NoiseLevel,TSlength,SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame)) %>% ungroup()
modelorder=unique(arrange(sims2_d, Classification, Model)$Model)
sims2_d$Model=factor(sims2_d$Model, levels=modelorder)
#set up results df
sims2_results=select(sims2_d, ID, SimNumber, Model)
#get hyperparameters (this takes a long time)
start=Sys.time()
sims2_results$hpar1=future_map(sims2_d$data, besthyper, y="Value", ylog=F, pgr="none")
end=Sys.time(); end-start
#get model output
sims2_results$modelresults1=map2(sims2_d$data, sims2_results$hpar1, smap_model_options, y="Value", model=1)
#get R2
sims2_d$R1a=map_dbl(sims2_results$modelresults1, ~.x$modelstats$R2abund)

sims2_d_hold=sims2_d
sims2_results_hold=sims2_results

save(sims2_d_hold, sims2_results_hold, file = "./data/sims2_results_hold.Rdata")

#50-100
sims2_d=select(sims, ID, Model, TimeStep, NoiseLevel, Classification, Sim.51:Sim.100) %>% 
  gather(SimNumber, Value, Sim.51:Sim.100) %>% 
  group_by(ID,Model,Classification,NoiseLevel,SimNumber) %>%  mutate(TSlength=length(Value)) %>% ungroup() %>% 
  group_by(ID,Model,Classification,NoiseLevel,TSlength,SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame)) %>% ungroup()
modelorder=unique(arrange(sims2_d, Classification, Model)$Model)
sims2_d$Model=factor(sims2_d$Model, levels=modelorder)
#set up results df
sims2_results=select(sims2_d, ID, SimNumber, Model)
#get hyperparameters (this take a long time)
start=Sys.time()
sims2_results$hpar1=future_map(sims2_d$data, besthyper, y="Value", ylog=F, pgr="none")
end=Sys.time(); end-start
#get model output
sims2_results$modelresults1=map2(sims2_d$data, sims2_results$hpar1, smap_model_options, y="Value", model=1)
#get R2
sims2_d$R1a=map_dbl(sims2_results$modelresults1, ~.x$modelstats$R2abund)

sims2_d_hold=rbind(sims2_d_hold, sims2_d)
sims2_results_hold=rbind(sims2_results_hold, sims2_results)

save(sims2_d_hold, sims2_results_hold, file = "./data/sims2_results_hold.Rdata")

sims2_d=sims2_d_hold
sims2_results=sims2_results_hold

#fit alt models for strictly positive data
# #1 and 3 are same, 2 and 4 are same
logmodels=c("PredatorPreyPeriodic", "PredatorPreyChaotic", "HostParParPeriodic", "HostParParChaotic")
sims2_log=filter(sims2_d, Model %in% logmodels)
sims2_log_results=filter(sims2_results, Model %in% logmodels)
sims2_log_results$hpar4=future_map(sims2_log$data, besthyper, y="Value", ylog=F, pgr="gr")
sims2_log_results$hpar5=future_map(sims2_log$data, besthyper, y="Value", ylog=T, pgr="gr")
sims2_log_results$modelresults4=map2(sims2_log$data, sims2_log_results$hpar4, smap_model_options, y="Value", model=4)
sims2_log_results$modelresults5=map2(sims2_log$data, sims2_log_results$hpar5, smap_model_options, y="Value", model=5)
sims2_log$R1a=map_dbl(sims2_log_results$modelresults1, ~.x$modelstats$R2abund)
sims2_log$R4a=map_dbl(sims2_log_results$modelresults4, ~.x$modelstats$R2abund)
sims2_log$R5a=map_dbl(sims2_log_results$modelresults5, ~.x$modelstats$R2abund)
sims2_log$bestR2=select(sims2_log,R1a,R4a,R5a) %>% apply(1,max)
sims2_log$bestmodel=select(sims2_log,R1a,R4a,R5a) %>% apply(1,which.max)
sims2_log_results$bestmodel=select(sims2_log,R1a,R4a,R5a) %>% apply(1,which.max)
sims2_log_results$modelresultsbest=cbind(select(sims2_log_results, modelresults1, modelresults4, modelresults5),sims2_log$bestmodel) %>% apply(1, function(x) {m=as.numeric(x["sims2_log$bestmodel"]); x[m][[1]]})

#join log models with others
sims2_results=left_join(sims2_results, select(sims2_log_results, ID, SimNumber, modelresultsbest, bestmodel), by=c("ID", "SimNumber"))
sims2_results$modelresultsbest=ifelse(sims2_d$Model %in% logmodels,
                                     sims2_results$modelresultsbest, 
                                     sims2_results$modelresults1)

#get stability
sims2_results$jacobians=map(sims2_results$modelresultsbest, getJacobians)
sims2_results$stability=map2(sims2_results$modelresultsbest, sims2_results$jacobians, getStability)
sims2_results$LEshift=map2(sims2_results$modelresultsbest, sims2_results$jacobians, LEshift)

#pull results
sims2_d$Ebest=map_dbl(sims2_results$modelresultsbest, ~.x$modelstats$E)
sims2_d$taubest=map_dbl(sims2_results$modelresultsbest, ~.x$modelstats$tau)
sims2_d$thetabest=map_dbl(sims2_results$modelresultsbest, ~.x$modelstats$theta)
sims2_d$R2best=map_dbl(sims2_results$modelresultsbest, ~.x$modelstats$R2abund)
sims2_d$gle=map_dbl(sims2_results$stability, ~.x$gle)
sims2_d$minmean=map_dbl(sims2_results$LEshift, ~.x$minmean)
sims2_d$minci=map_dbl(sims2_results$LEshift, ~.x$minci)

#join with first 20 reps
sims_d=rbind(sims_d, sims2_d)
sims_results=rbind(sims_results, sims2_results)

#regression method ####
sims_results$regLE=map2(sims_d$data, sims_results$modelresultsbest, regLE, y="Value")
sims_d$LEreg=map_dbl(sims_results$regLE, ~.x$LEreg)
sims_d$LEreg_se=map_dbl(sims_results$regLE, ~.x$LEreg_se)

#save results
#save(sims_d, sims_results, file = "./data/sims_results_update2.Rdata")
#save(sims_log, sims_log_results, sims2_log, sims2_log_results, file = "./data/sims_results_log_update.Rdata")

#export E and tau for other analyses
Eexport=spread(select(sims_d, ID:SimNumber, Ebest), SimNumber, Ebest) %>% 
  select(ID:TSlength, paste0("Sim.",1:100))
tauexport=spread(select(sims_d, ID:SimNumber, taubest), SimNumber, taubest) %>% 
  select(ID:TSlength, paste0("Sim.",1:100))
write.csv(Eexport,"./data/sims_test_E.csv", row.names = F)
write.csv(tauexport,"./data/sims_test_tau.csv", row.names = F)

#export results
sims_d$modelform=map_chr(sims_results$modelresultsbest, ~.x$form)
dexport=select(sims_d, ID:SimNumber, E=Ebest, tau=taubest, theta=thetabest, R2=R2best, modelform, LEmean=minmean, LEmin=minci, LEreg, LEreg_se)
write.csv(dexport,"./data/sims_test_results.csv", row.names = F)

#validation data ####
sims=read.csv("./data/ChaosMetaAnalysisSimulatedDataVALIDATION.csv")

sims_v=select(sims, ID, Model, TimeStep, NoiseLevel, TSlength=TimeSeriesLength, Classification, Sim.1:Sim.100) %>% 
  gather(SimNumber, Value, Sim.1:Sim.100) %>% 
  group_by(ID,Model,Classification,NoiseLevel,TSlength,SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame)) %>% ungroup()
modelorder=unique(arrange(sims_v, Classification, Model)$Model)
sims_v$Model=factor(sims_v$Model, levels=modelorder)

#set up results df
sims_vresults=select(sims_v, ID, SimNumber, Model)
#get hyperparameters (this takes a long time)
start=Sys.time()
sims_vresults$hpar1=future_map(sims_v$data, besthyper, y="Value", ylog=F, pgr="none")
end=Sys.time(); end-start
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
sims_vresults$stability=map2(sims_vresults$modelresultsbest, sims_vresults$jacobians, getStability)
sims_vresults$LEshift=map2(sims_vresults$modelresultsbest, sims_vresults$jacobians, LEshift)

#pull results
sims_v$Ebest=map_dbl(sims_vresults$modelresultsbest, ~.x$modelstats$E)
sims_v$taubest=map_dbl(sims_vresults$modelresultsbest, ~.x$modelstats$tau)
sims_v$thetabest=map_dbl(sims_vresults$modelresultsbest, ~.x$modelstats$theta)
sims_v$R2best=map_dbl(sims_vresults$modelresultsbest, ~.x$modelstats$R2abund)
sims_v$gle=map_dbl(sims_vresults$stability, ~.x$gle)
sims_v$minmean=map_dbl(sims_vresults$LEshift, ~.x$minmean)
sims_v$minci=map_dbl(sims_vresults$LEshift, ~.x$minci)

#regression method
sims_vresults$regLE=map2(sims_v$data, sims_vresults$modelresultsbest, regLE, y="Value")
sims_v$LEreg=map_dbl(sims_vresults$regLE, ~.x$LEreg)
sims_v$LEreg_se=map_dbl(sims_vresults$regLE, ~.x$LEreg_se)

#save results
#save(sims_v, sims_vresults, file = "./data/sims_validation2.Rdata")
#save(sims_vlog, sims_vlog_results, file = "./data/sims_validation_log.Rdata")

#export E and tau for other analyses
Eexport=spread(select(sims_v, ID:SimNumber, Ebest), SimNumber, Ebest) %>% 
  select(ID:TSlength, paste0("Sim.",1:100))
tauexport=spread(select(sims_v, ID:SimNumber, taubest), SimNumber, taubest) %>% 
  select(ID:TSlength, paste0("Sim.",1:100))
write.csv(Eexport,"./data/sims_validation_E.csv", row.names = F)
write.csv(tauexport,"./data/sims_validation_tau.csv", row.names = F)

#export results
sims_v$modelform=map_chr(sims_vresults$modelresultsbest, ~.x$form)
vexport=select(sims_v, ID:SimNumber, E=Ebest, tau=taubest, theta=thetabest, R2=R2best, modelform, LEmean=minmean, LEmin=minci, LEreg, LEreg_se)
write.csv(vexport,"./data/sims_validation_results.csv", row.names = F)

#more validation data ####
sims=read.csv("./data/SeasonalAR_data.csv")

sims_v2=select(sims, ID, Model, TimeStep, NoiseLevel, TSlength=TimeSeriesLength, Classification, Sim.1:Sim.100) %>% 
  gather(SimNumber, Value, Sim.1:Sim.100) %>% 
  group_by(ID,Model,Classification,NoiseLevel,TSlength,SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame)) %>% ungroup()
modelorder=unique(arrange(sims_v2, Classification, Model)$Model)
sims_v2$Model=factor(sims_v2$Model, levels=modelorder)

#set up results df
sims_vresults2=select(sims_v2, ID, SimNumber, Model)
#get hyperparameters (this takes a long time)
start=Sys.time()
sims_vresults2$hpar1=future_map(sims_v2$data, besthyper, y="Value", ylog=F, pgr="none")
end=Sys.time(); end-start
#get model output
sims_vresults2$modelresults1=map2(sims_v2$data, sims_vresults2$hpar1, smap_model_options, y="Value", model=1)
#get R2
sims_v2$R1a=map_dbl(sims_vresults2$modelresults1, ~.x$modelstats$R2abund)

#only 1 model
sims_vresults2$modelresultsbest=sims_vresults2$modelresults1

#get stability
sims_vresults2$jacobians=map(sims_vresults2$modelresultsbest, getJacobians)
sims_vresults2$stability=map2(sims_vresults2$modelresultsbest, sims_vresults2$jacobians, getStability)
sims_vresults2$LEshift=map2(sims_vresults2$modelresultsbest, sims_vresults2$jacobians, LEshift)

#pull results
sims_v2$Ebest=map_dbl(sims_vresults2$modelresultsbest, ~.x$modelstats$E)
sims_v2$taubest=map_dbl(sims_vresults2$modelresultsbest, ~.x$modelstats$tau)
sims_v2$thetabest=map_dbl(sims_vresults2$modelresultsbest, ~.x$modelstats$theta)
sims_v2$R2best=map_dbl(sims_vresults2$modelresultsbest, ~.x$modelstats$R2abund)
sims_v2$gle=map_dbl(sims_vresults2$stability, ~.x$gle)
sims_v2$minmean=map_dbl(sims_vresults2$LEshift, ~.x$minmean)
sims_v2$minci=map_dbl(sims_vresults2$LEshift, ~.x$minci)

#regression method
sims_vresults2$regLE=map2(sims_v2$data, sims_vresults2$modelresultsbest, regLE, y="Value")
sims_v2$LEreg=map_dbl(sims_vresults2$regLE, ~.x$LEreg)
sims_v2$LEreg_se=map_dbl(sims_vresults2$regLE, ~.x$LEreg_se)

save(sims_v2, sims_vresults2, file = "./data/sims_validation_forcedAR.Rdata")

#export E and tau for other analyses
Eexport=spread(select(sims_v2, ID:SimNumber, Ebest), SimNumber, Ebest) %>% 
  select(ID:TSlength, paste0("Sim.",1:100))
tauexport=spread(select(sims_v2, ID:SimNumber, taubest), SimNumber, taubest) %>% 
  select(ID:TSlength, paste0("Sim.",1:100))
write.csv(Eexport,"./data/sims_validation_E_forcedAR.csv", row.names = F)
write.csv(tauexport,"./data/sims_validation_tau_forcedAR.csv", row.names = F)

#export results
sims_v2$modelform=map_chr(sims_vresults2$modelresultsbest, ~.x$form)
vexport=select(sims_v2, ID:SimNumber, E=Ebest, tau=taubest, theta=thetabest, R2=R2best, modelform, LEmean=minmean, LEmin=minci, LEreg, LEreg_se)
write.csv(vexport,"./data/sims_validation_results_forcedAR.csv", row.names = F)


#merge with results from other methods ####

#test data 
#load("./data/sims_results_update.Rdata")
sims_d=read.csv("./data/sims_test_results.csv", stringsAsFactors = F)
modelorder=unique(arrange(sims_d, Classification, Model)$Model)
sims_d$Model=factor(sims_d$Model, levels=modelorder)
#reclass noise level for stochastic ts
sims_d$NoiseLevel2=ifelse(sims_d$NoiseLevel==0, 0.01, sims_d$NoiseLevel)
sims_d$Classification2=ifelse(sims_d$Classification=="chaotic", "chaotic", "not chaotic")
#class LEs regression method
sims_d$LEregmin=sims_d$LEreg-1.96*sims_d$LEreg_se
sims_d$LEregclass=ifelse(sims_d$LEregmin>0.01, "chaotic", "not chaotic")
#class LEs Jacobian method
sims_d$LEclass=ifelse(sims_d$LEmin>0.01, "chaotic", "not chaotic")
#class other methods
RQAclass=read.csv("./data/RQAclassification.csv") %>% gather(SimNumber, RQAclass, Sim.1:Sim.100)
PEclass=read.csv("./data/PEclassification.csv") %>% gather(SimNumber, PEclass, Sim.1:Sim.100)
HVAclass=read.csv("./data/HVAclassification.csv") %>% gather(SimNumber, HVAclass, Sim.1:Sim.100)
DTclass=read.csv("./data/DTclassification.csv") %>% gather(SimNumber, DTclass, Sim.1:Sim.100)
#join to main table
sims_d=left_join(sims_d, RQAclass) %>% left_join(PEclass) %>% 
  left_join(HVAclass) %>% left_join(DTclass)

#validation data
sims_v=read.csv("./data/sims_validation_results.csv", stringsAsFactors = F)
modelorderv=unique(arrange(sims_v, Classification, Model)$Model)
sims_v$Model=factor(sims_v$Model, levels=modelorderv)
#reclass noise level for stochastic ts
sims_v$NoiseLevel2=ifelse(sims_v$NoiseLevel==0, 0.01, sims_v$NoiseLevel)
sims_v$Classification2=ifelse(sims_v$Classification=="chaotic", "chaotic", "not chaotic")
#class LEs regression method
sims_v$LEregmin=sims_v$LEreg-1.96*sims_v$LEreg_se
sims_v$LEregclass=ifelse(sims_v$LEregmin>0.01, "chaotic", "not chaotic")
#class LEs jacobian method
sims_v$LEclass=ifelse(sims_v$LEmin>0.01, "chaotic", "not chaotic")
#class other methods
RQAclassv=read.csv("./data/RQAclassification_validation.csv") %>% gather(SimNumber, RQAclass, Sim.1:Sim.100)
PEclassv=read.csv("./data/PEclassification_validation.csv") %>% gather(SimNumber, PEclass, Sim.1:Sim.100)
HVAclassv=read.csv("./data/HVAclassification_validation.csv") %>% gather(SimNumber, HVAclass, Sim.1:Sim.100)
DTclassv=read.csv("./data/DTclassification_validation.csv") %>% gather(SimNumber, DTclass, Sim.1:Sim.100)
#join to main table
sims_v=left_join(sims_v, RQAclassv) %>% left_join(PEclassv) %>%
  left_join(HVAclassv) %>% left_join(DTclassv)

#more validation data
sims_v2=read.csv("./data/sims_validation_results_forcedAR.csv", stringsAsFactors = F)
modelorderv=unique(arrange(sims_v2, Classification, Model)$Model)
sims_v2$Model=factor(sims_v2$Model, levels=modelorderv)
#reclass noise level for stochastic ts
sims_v2$NoiseLevel2=ifelse(sims_v2$NoiseLevel==0, 0.01, sims_v2$NoiseLevel)
sims_v2$Classification2=ifelse(sims_v2$Classification=="chaotic", "chaotic", "not chaotic")
#class LEs regression method
sims_v2$LEregmin=sims_v2$LEreg-1.96*sims_v2$LEreg_se
sims_v2$LEregclass=ifelse(sims_v2$LEregmin>0.01, "chaotic", "not chaotic")
#class LEs jacobian method
sims_v2$LEclass=ifelse(sims_v2$LEmin>0.01, "chaotic", "not chaotic")
#class other methods
RQAclassv=read.csv("./data/RQAclassification_forcedAR.csv", stringsAsFactors = F) %>% gather(SimNumber, RQAclass, Sim.1:Sim.100)
PEclassv=read.csv("./data/PEclassification_forcedAR.csv", stringsAsFactors = F) %>% gather(SimNumber, PEclass, Sim.1:Sim.100)
HVAclassv=read.csv("./data/HVAclassification_forcedAR.csv", stringsAsFactors = F) %>% gather(SimNumber, HVAclass, Sim.1:Sim.100)
DTclassv=read.csv("./data/DTclassification_forcedAR.csv", stringsAsFactors = F) %>% gather(SimNumber, DTclass, Sim.1:Sim.100)
#join to main table
sims_v2=left_join(sims_v2, RQAclassv) %>% left_join(PEclassv) %>%
  left_join(HVAclassv) %>% left_join(DTclassv)

#write data
write.csv(sims_d, "./data/sims_test_results_allmethods.csv", row.names = F)
write.csv(sims_v, "./data/sims_validation_results_allmethods.csv", row.names = F)
write.csv(sims_v2, "./data/sims_validation_results_forcedAR_allmethods.csv", row.names = F)

#test with known LE ####
sims=read.csv("./data/ModelsWithKnownLEs.csv")
sims_d=sims %>%  
  gather(SimNumber, Value, Sim.1:Sim.4) %>% 
  group_by(SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame))
sims_d$modelresults1=map(sims_d$data, smap_model_options, y="Value", model=1)
sims_d$jacobians=map(sims_d$modelresults1, getJacobians)
sims_d$LEshift=map2(sims_d$modelresults1, sims_d$jacobians, LEshift)
sims_d$minci=map_dbl(sims_d$LEshift, ~.x$minci)
sims_d$R2best=map_dbl(sims_d$modelresults1, ~.x$modelstats$R2abund)

sims_d$regLE=map2(sims_d$data, sims_d$modelresults1, regLE, y="Value")
sims_d$LEreg=map_dbl(sims_d$regLE, ~.x$LEreg)
sims_d$LEreg_se=map_dbl(sims_d$regLE, ~.x$LEreg_se)
sims_d$LEregmin=sims_d$LEreg-1.96*sims_d$LEreg_se

