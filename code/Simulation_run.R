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
#save(sims_d, sims_results, file = "./data/sims_results_update.Rdata")
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
#save(sims_v, sims_vresults, file = "./data/sims_validation.Rdata")
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

#testing ####
ser_or=sims_test[4,]$data[[1]]$Value

# #best hyperparameters
# sims_results$hyperpars=map(sims_d$data, besthyper, y="Value")
# sims_d=cbind(sims_d, bind_rows(sims_results$hyperpars))
# sims_d2=rbind(sims_d2, sims_d)

sims_test=filter(sims_d, Model %in% c("HostParParPeriodic") & TSlength==25)
sims_test=filter(sims_d, Model %in% c("HostParParChaotic"))
sims_test_results=filter(sims_results, ID %in% sims_test$ID)

sims_test_results$hpar4=map(sims_test$data, besthyper, y="Value", ylog=F, pgr="gr",returntable=T)
sims_test_results$hpar1=map(sims_test$data, besthyper, y="Value", ylog=F, pgr="none",returntable=T)

res1=cbind(select(sims_test, ID:SimNumber), map_df(sims_test_results$hpar1, function(x) {
  model_selection=arrange(x, tau, theta, E)
  #model_selection=arrange(x, E, tau, theta)
  bestE=model_selection$E[which.min(model_selection$error)]
  bestTau=model_selection$tau[which.min(model_selection$error)]
  bestTheta=model_selection$theta[which.min(model_selection$error)]
  return(data.frame(bestE=bestE, bestTau=bestTau, bestTheta=bestTheta))
}))

ggplot(res1, aes(x=bestE, y=bestTau, color=Model)) +
  facet_grid(TSlength~NoiseLevel) + 
  geom_point(alpha=0.5, size=3) + theme_bw()

sims_test_results$hparsum4=map(sims_test_results$hpar4, function(x) {
  model_selection=arrange(x, tau, theta, E) #%>% filter(tau==1)
  #model_selection=arrange(x, E, tau, theta)
  bestE=model_selection$E[which.min(model_selection$error)]
  bestTau=model_selection$tau[which.min(model_selection$error)]
  bestTheta=model_selection$theta[which.min(model_selection$error)]
  return(data.frame(bestE=bestE, bestTau=bestTau, bestTheta=bestTheta, Emax=max(model_selection$E), taumax=max(model_selection$tau)))
})
sims_test_results$modelresults4=map2(sims_test$data, sims_test_results$hparsum4, smap_model_options, y="Value", model=4)
sims_test_results$jacobians4=map(sims_test_results$modelresults4, getJacobians)
sims_test_results$stability4=map2(sims_test_results$modelresults4, sims_test_results$jacobians4, getStability)
sims_test$gle4=map_dbl(sims_test_results$stability4, ~.x$gle)
sims_test$Ebest4=map_dbl(sims_test_results$modelresults4, ~.x$modelstats$E)
sims_test$taubest4=map_dbl(sims_test_results$modelresults4, ~.x$modelstats$tau)
sims_test$thetabest4=map_dbl(sims_test_results$modelresults4, ~.x$modelstats$theta)
sims_test$R4a=map_dbl(sims_test_results$modelresults4, ~.x$modelstats$R2abund)
sims_test$gle_class4=ifelse(sims_test$gle4>0.05, "chaotic", "not chaotic")
sims_test$LEshift4=map2(sims_test_results$modelresults4, sims_test_results$jacobians4, LEshift)
sims_test$minmean=map_dbl(sims_test$LEshift4, ~.x$minmean)
sims_test$minci=map_dbl(sims_test$LEshift4, ~.x$minci)
sims_test$LEs_class4=ifelse(sims_test$minci>0.05, "chaotic", "not chaotic")

sims_test_results$hparsum1=map(sims_test_results$hpar1, function(x) {
  model_selection=arrange(x, tau, theta, E) #%>% filter(tau==1)
  #model_selection=arrange(x, E, tau, theta)
  bestE=model_selection$E[which.min(model_selection$error)]
  bestTau=model_selection$tau[which.min(model_selection$error)]
  bestTheta=model_selection$theta[which.min(model_selection$error)]
  return(data.frame(bestE=bestE, bestTau=bestTau, bestTheta=bestTheta, Emax=max(model_selection$E), taumax=max(model_selection$tau)))
})
sims_test_results$modelresults1=map2(sims_test$data, sims_test_results$hparsum1, smap_model_options, y="Value", model=1)
sims_test_results$jacobians1=map(sims_test_results$modelresults1, getJacobians)
sims_test_results$stability1=map2(sims_test_results$modelresults1, sims_test_results$jacobians1, getStability)
sims_test$gle1=map_dbl(sims_test_results$stability1, ~.x$gle)
sims_test$Ebest1=map_dbl(sims_test_results$modelresults1, ~.x$modelstats$E)
sims_test$taubest1=map_dbl(sims_test_results$modelresults1, ~.x$modelstats$tau)
sims_test$thetabest1=map_dbl(sims_test_results$modelresults1, ~.x$modelstats$theta)
sims_test$R1a=map_dbl(sims_test_results$modelresults1, ~.x$modelstats$R2abund)
sims_test$gle_class1=ifelse(sims_test$gle1>0.05, "chaotic", "not chaotic")
sims_test$LEshift1=map2(sims_test_results$modelresults1, sims_test_results$jacobians1, LEshift)
sims_test$minmean1=map_dbl(sims_test$LEshift1, ~.x$minmean)
sims_test$minci1=map_dbl(sims_test$LEshift1, ~.x$minci)
sims_test$LEs_class1=ifelse(sims_test$minci1>0.05, "chaotic", "not chaotic")

sims_summary=sims_test %>% select(-data) %>% 
  group_by(Classification,NoiseLevel,TSlength, Model) %>% 
  summarize(gle4_pp=length(which(gle4>0.0))/length(gle4),
            gle4_pp.01=length(which(gle4>0.01))/length(gle4),
            gle4_pp.05=length(which(gle4>0.05))/length(gle4),
            minci_pp=length(which(minci>0.0))/length(minci),
            minci_pp.01=length(which(minci>0.01))/length(minci),
            minci_pp.05=length(which(minci>0.05))/length(minci),
            minci1_pp=length(which(minci1>0.0))/length(minci1),
            minci1_pp.01=length(which(minci1>0.01))/length(minci1),
            minci1_pp.05=length(which(minci1>0.05))/length(minci1))
ggplot(sims_summary, aes(x=Model, y=minci1_pp.01, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel) + ggtitle("tau, theta, E") +
  geom_bar(stat = "identity") + theme_bw() + xlabvert

ggplot(sims_test, aes(x=Model, y=gle4, color=Classification)) +
  facet_grid(TSlength~NoiseLevel) + geom_hline(yintercept = 0) + geom_hline(yintercept = 0.05, lty=2) +
  geom_point(alpha=0.3, size=2) + theme_bw() + xlabvert + legalpha

#save(sims_test, sims_test_results, file = "./data/sims_results_problemmodel.Rdata")


sims_test_results$modelresults1=map(sims_test$data, smap_model_options, y="Value", model=1)
sims_test_results$modelresults4=map(sims_test$data, smap_model_options, y="Value", model=4)
sims_test_results$modelresults1b=pmap(list(sims_test$data, taufix=sims_test$Tau), smap_model_options, y="Value", model=1)
sims_test_results$modelresults1c=pmap(list(sims_test$data, Efix=sims_test$E, taufix=sims_test$Tau), smap_model_options, y="Value", model=1)
sims_test_results$modelresults1d=map(sims_test$data, smap_model_options, y="Value", model=1, taufix=2, Efix=3)
sims_test_results$jacobians1=map(sims_test_results$modelresults1, getJacobians)
sims_test_results$jacobians4=map(sims_test_results$modelresults4, getJacobians)
sims_test_results$jacobians1b=map(sims_test_results$modelresults1b, getJacobians)
sims_test_results$jacobians1c=map(sims_test_results$modelresults1c, getJacobians)
sims_test_results$jacobians1d=map(sims_test_results$modelresults1d, getJacobians)
sims_test_results$stability1=map2(sims_test_results$modelresults1, sims_test_results$jacobians1, getStability)
sims_test_results$stability4=map2(sims_test_results$modelresults4, sims_test_results$jacobians4, getStability)
sims_test_results$stability1b=map2(sims_test_results$modelresults1b, sims_test_results$jacobians1b, getStability)
sims_test_results$stability1c=map2(sims_test_results$modelresults1c, sims_test_results$jacobians1c, getStability)
sims_test_results$stability1d=map2(sims_test_results$modelresults1d, sims_test_results$jacobians1d, getStability)
sims_test$gle1=map_dbl(sims_test_results$stability1, ~.x$gle)
sims_test$gle4=map_dbl(sims_test_results$stability4, ~.x$gle)
sims_test$gle1b=map_dbl(sims_test_results$stability1b, ~.x$gle)
sims_test$gle1c=map_dbl(sims_test_results$stability1c, ~.x$gle)
sims_test$gle1d=map_dbl(sims_test_results$stability1d, ~.x$gle)
sims_test$Ebest1=map_dbl(sims_test_results$modelresults1, ~.x$modelstats$E)
sims_test$Ebest4=map_dbl(sims_test_results$modelresults4, ~.x$modelstats$E)
sims_test$Ebest1b=map_dbl(sims_test_results$modelresults1b, ~.x$modelstats$E)
sims_test$Ebest1c=map_dbl(sims_test_results$modelresults1c, ~.x$modelstats$E)
sims_test$Ebest1d=map_dbl(sims_test_results$modelresults1d, ~.x$modelstats$E)
sims_test$taubest1=map_dbl(sims_test_results$modelresults1, ~.x$modelstats$tau)
sims_test$taubest4=map_dbl(sims_test_results$modelresults4, ~.x$modelstats$tau)
sims_test$taubest1b=map_dbl(sims_test_results$modelresults1b, ~.x$modelstats$tau)
sims_test$thetabest1=map_dbl(sims_test_results$modelresults1, ~.x$modelstats$theta)
sims_test$thetabest4=map_dbl(sims_test_results$modelresults4, ~.x$modelstats$theta)
sims_test$thetabest1d=map_dbl(sims_test_results$modelresults1d, ~.x$modelstats$theta)
sims_test$R1a=map_dbl(sims_test_results$modelresults1, ~.x$modelstats$R2abund)
sims_test$R4a=map_dbl(sims_test_results$modelresults4, ~.x$modelstats$R2abund)
sims_test$gle_class4=ifelse(sims_test$gle4>0.01, "chaotic", "not chaotic")

sims_summary=sims_test %>% select(-data) %>% 
  group_by(Classification,NoiseLevel,TSlength, Model) %>% 
  summarize(gle4_pp=length(which(gle4>0.0))/length(gle4),
            gle4_pp.01=length(which(gle4>0.01))/length(gle4),
            gle4_pp.05=length(which(gle4>0.05))/length(gle4))
            # gle1b_pp=length(which(gle1b>0.0))/length(gle1b),
            # gle1b_pp.01=length(which(gle1b>0.01))/length(gle1b),
            # gle1c_pp=length(which(gle1c>0.0))/length(gle1c),
            # gle1c_pp.01=length(which(gle1c>0.01))/length(gle1c),
            # gle1d_pp=length(which(gle1d>0.0))/length(gle1d),
            # gle1d_pp.01=length(which(gle1d>0.01))/length(gle1d))
ggplot(sims_summary, aes(x=factor(NoiseLevel), y=gle4_pp.05)) +
  facet_grid(TSlength~.) + 
  geom_bar(stat = "identity") + theme_bw()
ggplot(sims_summary, aes(x=Model, y=gle4_pp.05, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel) + ggtitle("tau, E, theta") +
  geom_bar(stat = "identity") + theme_bw() + xlabvert


sims_test_results$modelresults2=map(sims_test$data, smap_model, y="Value", ylog=T)
sims_test_results$modelresults5=map(sims_test$data, smap_model, y="Value", ylog=F, taufix=2, Efix=4)
sims_test_results$jacobians1=map(sims_test_results$modelresults1, getJacobians)
sims_test_results$jacobians2=map(sims_test_results$modelresults2, getJacobians)
sims_test_results$jacobians5=map(sims_test_results$modelresults5, getJacobians)
sims_test_results$stability1=map2(sims_test_results$modelresults1, sims_test_results$jacobians1, getStability)
sims_test_results$stability2=map2(sims_test_results$modelresults2, sims_test_results$jacobians2, getStability)
sims_test_results$stability5=map2(sims_test_results$modelresults5, sims_test_results$jacobians5, getStability)

sims_test$LEshift1=map2_dbl(sims_test_results$modelresults1, sims_test_results$jacobians1, LEshift)
sims_test$LEshift2=map2_dbl(sims_test_results$modelresults2, sims_test_results$jacobians2, LEshift)
sims_test$LEshift5=map2_dbl(sims_test_results$modelresults5, sims_test_results$jacobians5, LEshift)
sims_test$LEshiftbest=map2_dbl(sims_test_results$modelresults2, sims_test_results$jacobians2, LEshift)
sims_test$gle1=map_dbl(sims_test_results$stability1, ~.x$gle)
sims_test$gle2=map_dbl(sims_test_results$stability2, ~.x$gle)
sims_test$gle5=map_dbl(sims_test_results$stability5, ~.x$gle)
sims_test$Ebest1=map_dbl(sims_test_results$modelresults1, ~.x$modelstats$E)
sims_test$Ebest2=map_dbl(sims_test_results$modelresults2, ~.x$modelstats$E)
sims_test$Ebest5=map_dbl(sims_test_results$modelresults5, ~.x$modelstats$E)
sims_test$taubest1=map_dbl(sims_test_results$modelresults1, ~.x$modelstats$tau)
sims_test$taubest2=map_dbl(sims_test_results$modelresults2, ~.x$modelstats$tau)
sims_test$thetabest1=map_dbl(sims_test_results$modelresults1, ~.x$modelstats$theta)
sims_test$thetabest2=map_dbl(sims_test_results$modelresults2, ~.x$modelstats$theta)
sims_test$R1a=map_dbl(sims_test_results$modelresults1, ~.x$modelstats$R2abund)
sims_test$R2a=map_dbl(sims_test_results$modelresults2, ~.x$modelstats$R2abund)
sims_test$R5a=map_dbl(sims_test_results$modelresults5, ~.x$modelstats$R2abund)

sims_test$bestmodel=select(sims_test,R1a,R2a) %>% apply(1,which.max)
sims_test$R2best=select(sims_test,R1a,R2a) %>% apply(1,max)
sims_test$glebest=apply(sims_test, 1, FUN=function(x) {unlist(x[paste0("gle",x$bestmodel)])})
sims_test$LEshiftbest=apply(sims_test, 1, FUN=function(x) {unlist(x[paste0("LEshift",x$bestmodel)])})
sims_test$Ebest=apply(sims_test, 1, FUN=function(x) {unlist(x[paste0("Ebest",x$bestmodel)])})
sims_test$taubest=apply(sims_test, 1, FUN=function(x) {unlist(x[paste0("taubest",x$bestmodel)])})

sims_summary=sims_test %>% select(-data) %>% 
  group_by(Classification,NoiseLevel,TSlength, Model) %>% 
  summarize(glebest_pp=length(which(glebest>0.0))/length(gle1),
            glebest_pp.01=length(which(glebest>0.01))/length(gle1),
            glebest_pp.05=length(which(glebest>0.05))/length(gle1),
            LEshiftbest_pp=length(which(LEshiftbest>0))/length(LEshift1),
            LEshiftbest_pp.01=length(which(LEshiftbest>0.01))/length(LEshift1),
            LEshiftbest_pp.05=length(which(LEshiftbest>0.05))/length(LEshift1),
            gle2_pp=length(which(gle2>0.0))/length(gle2),
            gle2_pp.01=length(which(gle2>0.01))/length(gle2),
            gle2_pp.05=length(which(gle2>0.05))/length(gle2),
            LEshift2_pp=length(which(LEshift2>0))/length(LEshift2),
            LEshift2_pp.01=length(which(LEshift2>0.01))/length(LEshift2),
            LEshift2_pp.05=length(which(LEshift2>0.05))/length(LEshift2))


ggplot(sims_test, aes(x=Model, y=taubest, color=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point(alpha=0.3) + theme_bw() + xlabvert + legalpha

ggplot(sims_summary, aes(x=Model, y=gle1_pp.01, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + 
  geom_bar(stat = "identity") + theme_bw() + xlabvert
ggplot(sims_summary, aes(x=factor(NoiseLevel), y=gle1d_pp.01)) +
  facet_grid(TSlength~., scales = "free_y") + 
  geom_bar(stat = "identity") + theme_bw()
#

sims_zero=filter(sims_test, Model=="HostParParPeriodic" & NoiseLevel==0 & TSlength>25)
sims_zero_results=filter(sims_test_results, ID %in% sims_zero$ID)

sims_zero_results$modelresults5=map(sims_zero$data, smap_model, y="Value", ylog=F, taufix=1, Efix=3)
sims_zero_results$jacobians5=map(sims_zero_results$modelresults5, getJacobians)
sims_zero_results$stability5=map2(sims_zero_results$modelresults5, sims_zero_results$jacobians5, getStability)
sims_zero$gle5=map_dbl(sims_zero_results$stability5, ~.x$gle)
sims_zero$LEshift5=map2_dbl(sims_zero_results$modelresults5, sims_zero_results$jacobians5, LEshift)

ggplot(filter(sims_test, Model=="HostParParPeriodic"), aes(x=Ebest, y=taubest)) +
  facet_grid(TSlength~NoiseLevel) +
  geom_point(x=4, y=2, color="red", size=2)+
  geom_point(alpha=0.3) + theme_bw()

# rescale=function(Model, data) {
#   if(Model %in% c("PredatorPreyPeriodic", "PredatorPreyChaotic", "HostParParPeriodic", "HostParParChaotic")) {
#     data$Value=log(data$Value)
#   }
#   return(as.data.frame(data))
# }
# sims_test$data_rescale=map2(sims_test$Model, sims_test$data, rescale)

sims_test_results$jacobians1=map(sims_test_results$modelresults2, getJacobians)
sims_test_results$stability1=map2(sims_test_results$modelresults2, sims_test_results$jacobians1, getStability)

sims_test_results$modelresults1b=map_if(sims_test$data_rescale, sims_test$Model %in% c("PredatorPreyPeriodic", "PredatorPreyChaotic", "HostParParPeriodic", "HostParParChaotic"),
                                        smap_model, y="Value", ylog=F)
sims_test_results$modelresults1b=ifelse(sims_test$Model %in% c("PredatorPreyPeriodic", "PredatorPreyChaotic", "HostParParPeriodic", "HostParParChaotic"),
                                        sims_test_results$modelresults1b, sims_test_results$modelresults1)
sims_test_results$jacobians1b=map(sims_test_results$modelresults1b, getJacobians)
sims_test_results$stability1b=map2(sims_test_results$modelresults1b, sims_test_results$jacobians1b, getStability)
sims_test$LEshift=map2_dbl(sims_test_results$modelresults1b, sims_test_results$jacobians1b, LEshift)
sims_test$gle=map_dbl(sims_test_results$stability1b, ~.x$gle)
sims_test$R1a=map_dbl(sims_test_results$modelresults1, ~.x$modelstats$R2abund)

# #external E and tau
# simE=read.csv("./data/ESimulatedNEW.csv")
# simtau=read.csv("./data/TauSimulatedNEW.csv")
# sims_E=select(simE, ID, Sim.1:Sim.5) %>% 
#   gather(SimNumber, E, Sim.1:Sim.5)
# sims_tau=select(simtau, ID, Sim.1:Sim.5) %>% 
#   gather(SimNumber, Tau, Sim.1:Sim.5)
# sims_Etau=left_join(sims_E, sims_tau, by=c("ID", "SimNumber"))
