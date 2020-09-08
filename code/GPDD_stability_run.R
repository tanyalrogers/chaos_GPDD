# GPDD stability analysis - run

library("rgpdd")
library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")
library("purrr")

source("./code/GPDD_stability_functions.R")
load(file = "./data/gpdd_d.Rdata")

# focal=gpdd_d$MainID
# gpdd_results=filter(gpdd_results, MainID %in% focal)
# gpdd_d=filter(gpdd_d, MainID %in% focal)

# fixed E and tau values to use
# Etau=read.csv("./data/EandTauEmpiricalNEW.csv")
# #gpdd_d=left_join(gpdd_d, Etau, by="MainID")
# gpdd_d=cbind(gpdd_d, Etau[,2:3])

gpdd_results=select(gpdd_d, MainID)

# #hyperpars
# gpdd_results$hyperpars=map(gpdd_d$data_rescale, besthyper, y="PopRescale")
# gpdd_results$hyperparslog=map(gpdd_d$data_rescale, besthyper, y="PopRescale_log")
# gpdd_d=cbind(gpdd_d, bind_rows(gpdd_results$hyperparslog))

#fit models (unconstrained E, tau)
#ut-ut
gpdd_results$modelresults1=map(gpdd_d$data_rescale, smap_model_options, y="PopRescale", model=1)
#log-log
gpdd_results$modelresults2=map(gpdd_d$data_rescale, smap_model_options, y="PopRescale", model=2)
#fd-ut
gpdd_results$modelresults3=map(gpdd_d$data_rescale, smap_model_options, y="PopRescale", model=3)
#gr-ut
gpdd_results$modelresults4=map(gpdd_d$data_rescale, smap_model_options, y="PopRescale", model=4)
#gr-log
gpdd_results$modelresults5=map(gpdd_d$data_rescale, smap_model_options, y="PopRescale", model=5)

#fit models (fixed E, tau)
# #ut-ut
# gpdd_results$modelresults1=pmap(list(data=gpdd_d$data_rescale, taufix=gpdd_d$Tau), smap_model_options, y="PopRescale", model=1)
# #log-log
# gpdd_results$modelresults2=pmap(list(data=gpdd_d$data_rescale, taufix=gpdd_d$Tau), smap_model_options, y="PopRescale", model=2)
# #fd-ut
# gpdd_results$modelresults3=pmap(list(data=gpdd_d$data_rescale, Efix=gpdd_d$E, taufix=gpdd_d$Tau), smap_model_options, y="PopRescale", model=3)
# #gr-ut
# gpdd_results$modelresults4=pmap(list(data=gpdd_d$data_rescale, Efix=gpdd_d$E, taufix=gpdd_d$Tau), smap_model_options, y="PopRescale", model=4)
# #gr-log
# gpdd_results$modelresults5=pmap(list(data=gpdd_d$data_rescale, Efix=gpdd_d$E, taufix=gpdd_d$Tau), smap_model_options, y="PopRescale", model=5)

#pull out R2 values for each model
gpdd_d$R1m=map_dbl(gpdd_results$modelresults1, ~.x$modelstats$R2model)
gpdd_d$R2m=map_dbl(gpdd_results$modelresults2, ~.x$modelstats$R2model)
gpdd_d$R3m=map_dbl(gpdd_results$modelresults3, ~.x$modelstats$R2model)
gpdd_d$R4m=map_dbl(gpdd_results$modelresults4, ~.x$modelstats$R2model)
gpdd_d$R5m=map_dbl(gpdd_results$modelresults5, ~.x$modelstats$R2model)
gpdd_d$R1a=map_dbl(gpdd_results$modelresults1, ~.x$modelstats$R2abund)
gpdd_d$R2a=map_dbl(gpdd_results$modelresults2, ~.x$modelstats$R2abund)
gpdd_d$R3a=map_dbl(gpdd_results$modelresults3, ~.x$modelstats$R2abund)
gpdd_d$R4a=map_dbl(gpdd_results$modelresults4, ~.x$modelstats$R2abund)
gpdd_d$R5a=map_dbl(gpdd_results$modelresults5, ~.x$modelstats$R2abund)
# gpdd_d=cbind(gpdd_d,map_df(gpdd_d$modelresults1, ~.x$modelstats))

#best model
gpdd_d$bestmodel=select(gpdd_d,R3a,R4a,R5a) %>% apply(1,which.max)
gpdd_d$bestR2=select(gpdd_d,R3a,R4a,R5a) %>% apply(1,max)
gpdd_d$bestR2m=select(gpdd_d,R3m,R4m,R5m) %>% apply(1,max)
gpdd_results$modelresultsbest=cbind(select(gpdd_results, modelresults3, modelresults4, modelresults5),gpdd_d$bestmodel) %>% apply(1, function(x) {m=as.numeric(x["gpdd_d$bestmodel"]); x[m][[1]]})
gpdd_d$E=map_dbl(gpdd_results$modelresultsbest, ~.x$modelstats$E)
gpdd_d$tau=map_dbl(gpdd_results$modelresultsbest, ~.x$modelstats$tau)
gpdd_d$theta=map_dbl(gpdd_results$modelresultsbest, ~.x$modelstats$theta)

#jacobians and stability of best model
gpdd_results$jacobians=map(gpdd_results$modelresultsbest, getJacobians)
gpdd_results$stability=map2(gpdd_results$modelresultsbest, gpdd_results$jacobians, getStability)
gpdd_results$LEshift=map2(gpdd_results$modelresultsbest, gpdd_results$jacobians, LEshift)
#LEs of best model
gpdd_d$gle=map_dbl(gpdd_results$stability, ~.x$gle)
gpdd_d$minci=map_dbl(gpdd_results$LEshift, ~.x$minci)
gpdd_d$lle_pp=map_dbl(gpdd_results$stability, ~.x$lle_pp)
#LE signs
gpdd_d$glesign=ifelse(gpdd_d$gle>0.01, "chaotic", "not chaotic")
gpdd_d$LEshiftsign=ifelse(gpdd_d$minci>0.01, "chaotic", "not chaotic")

length(which(gpdd_d$gle>0.01))/length(which(!is.na(gpdd_d$gle)))
length(which(gpdd_d$minci>0.01))/length(which(!is.na(gpdd_d$minci)))

#predictability
predthreshold=0.2
gpdd_d$predictable=ifelse(gpdd_d$bestR2>predthreshold, "yes","no")
gpdd_d$predictable_gr=ifelse(gpdd_d$bestR2m>predthreshold, "yes","no")
gpdd_d$predictable_ag=ifelse(gpdd_d$predictable=="yes" & gpdd_d$predictable_gr=="yes", "ag",
                              ifelse(gpdd_d$predictable=="yes" & gpdd_d$predictable_gr=="no", "a",
                                     ifelse(gpdd_d$predictable=="no" & gpdd_d$predictable_gr=="yes", "g", "none")))

# #verify 1&3 and 2&5 are the same
# #get jacobians
# gpdd_results2$jacobians1=map(gpdd_results$modelresults1, getJacobians)
# gpdd_results2$jacobians2=map(gpdd_results$modelresults2, getJacobians)
# gpdd_results2$jacobians3=map(gpdd_results$modelresults3, getJacobians)
# gpdd_results2$jacobians4=map(gpdd_results$modelresults4, getJacobians)
# gpdd_results2$jacobians5=map(gpdd_results$modelresults5, getJacobians)
# #get stability
# gpdd_results2$stability1=map2(gpdd_results$modelresults1, gpdd_results2$jacobians1, getStability)
# gpdd_results2$stability2=map2(gpdd_results$modelresults2, gpdd_results2$jacobians2, getStability)
# gpdd_results2$stability3=map2(gpdd_results$modelresults3, gpdd_results2$jacobians3, getStability)
# gpdd_results2$stability4=map2(gpdd_results$modelresults4, gpdd_results2$jacobians4, getStability)
# gpdd_results2$stability5=map2(gpdd_results$modelresults5, gpdd_results2$jacobians5, getStability)
# 
# gpdd_d2$gle1=map_dbl(gpdd_results2$stability1, ~.x$gle)
# gpdd_d2$gle2=map_dbl(gpdd_results2$stability2, ~.x$gle)
# gpdd_d2$gle3=map_dbl(gpdd_results2$stability3, ~.x$gle)
# gpdd_d2$gle4=map_dbl(gpdd_results2$stability4, ~.x$gle)
# gpdd_d2$gle5=map_dbl(gpdd_results2$stability5, ~.x$gle)

#convert to common timescale
gpdd_d$gle_mo=gpdd_d$gle/timescale_mo(gpdd_d$SamplingInterval, 1)
gpdd_d$gle_gen=gpdd_d$gle_mo*gpdd_d$MinAge_mo
gpdd_d$les_mo=gpdd_d$minci/timescale_mo(gpdd_d$SamplingInterval, 1)
gpdd_d$les_gen=gpdd_d$les_mo*gpdd_d$MinAge_mo

# #regression method
# gpdd_results$regLE=map(gpdd_d$data_rescale, regLE, y="PopRescale")
# gpdd_d$LEreg=map_dbl(gpdd_results$regLE, ~.x$LEreg)
# gpdd_d$LEreg_se=map_dbl(gpdd_results$regLE, ~.x$LEreg_se)

#fix E to 1 for best model
#tau unconstrained
gpdd_results$LE1d=map2(gpdd_d$data_rescale, gpdd_d$bestmodel+2, LE1d, y="PopRescale")
gpdd_results$jacobians1d=map(gpdd_results$LE1d, ~.x$jacobians)
gpdd_results$modelresults1d=map(gpdd_results$LE1d, ~.x$modelresults)
gpdd_results$LEshift1d=map2(gpdd_results$modelresults1d, gpdd_results$jacobians1d, LEshift)

# #tau fixed
# gpdd_results$LE1d=pmap(list(data=gpdd_d$data_rescale, bestmodel=gpdd_d$bestmodel, taufix=gpdd_d$Tau), LE1d, y="PopRescale")

gpdd_d$gle1d=map_dbl(gpdd_results$LE1d, ~.x$stability$gle)
gpdd_d$minci1d=map_dbl(gpdd_results$LEshift1d, ~.x$minci)
gpdd_d$LEshiftsign1d=ifelse(gpdd_d$minci1d>0.01, "chaotic", "not chaotic")

length(which(gpdd_d$gle1d>0.01))/length(which(!is.na(gpdd_d$gle1d)))
length(which(gpdd_d$minci1d>0.01))/length(which(!is.na(gpdd_d$minci1d)))

#prediction slope for best model



save(gpdd_d, gpdd_results, file = "./Data/gpdd_results_update.Rdata")

# #increase threshold
# predthreshold=0.2
# gpdd_d$predictable=select(gpdd_d,R1a,R2a,R3a,R4a,R5a) %>% apply(1,FUN = function(x) ifelse(any(x>predthreshold), "yes","no"))
# gpdd_d$predictable_gr=select(gpdd_d,R3m,R4m,R5m) %>% apply(1,FUN = function(x) ifelse(any(x>predthreshold), "yes","no"))
# gpdd_d$predictable_ag=ifelse(gpdd_d$predictable=="yes" & gpdd_d2$predictable_gr=="yes", "ag",
#                              ifelse(gpdd_d2$predictable=="yes" & gpdd_d2$predictable_gr=="no", "a",
#                                     ifelse(gpdd_d2$predictable=="no" & gpdd_d2$predictable_gr=="yes", "g", "none")))
# gpdd_d3=filter(gpdd_d, predictable_ag!="none")
# gpdd_results3=filter(gpdd_results, MainID %in% unique(gpdd_d3$MainID))
