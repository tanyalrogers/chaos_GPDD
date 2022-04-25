# Applies Jacobian LE chaos detection method to GPDD
# Includes selection of E and tau used by RQA

library(rEDM)
library(dplyr)
library(tidyr)
library(purrr)

source("./code/Methods/LE_ChaosDetectionMethods.R")

gpdd_d=read.csv("./data/gpdd_ts_metadata.csv")
gpdd_d_ts=read.csv("./data/gpdd_timeseries.csv")
gpdd_d_ts=gpdd_d_ts %>% group_by(MainID) %>% nest(.key="data_rescale") %>% 
  mutate(data_rescale=map(data_rescale, as.data.frame))
gpdd_d=left_join(gpdd_d, gpdd_d_ts, by="MainID")

gpdd_results=select(gpdd_d, MainID)

#### Jacobian LE Method ####

#fit models
#models 1 and 2 were dropped because the results of 1&3 and 2&5 are identical
gpdd_results$modelresults3=map(gpdd_d$data_rescale, smap_model_options, y="PopRescale", model=3) #fd-ut
gpdd_results$modelresults4=map(gpdd_d$data_rescale, smap_model_options, y="PopRescale", model=4) #gr-ut
gpdd_results$modelresults5=map(gpdd_d$data_rescale, smap_model_options, y="PopRescale", model=5) #gr-log

#pull out R2 values for each model
gpdd_d$R3m=map_dbl(gpdd_results$modelresults3, ~.x$modelstats$R2model)
gpdd_d$R4m=map_dbl(gpdd_results$modelresults4, ~.x$modelstats$R2model)
gpdd_d$R5m=map_dbl(gpdd_results$modelresults5, ~.x$modelstats$R2model)
gpdd_d$R3a=map_dbl(gpdd_results$modelresults3, ~.x$modelstats$R2abund)
gpdd_d$R4a=map_dbl(gpdd_results$modelresults4, ~.x$modelstats$R2abund)
gpdd_d$R5a=map_dbl(gpdd_results$modelresults5, ~.x$modelstats$R2abund)

#get best model of the 3 forms
gpdd_d$bestmodel=select(gpdd_d,R3a,R4a,R5a) %>% apply(1,which.max)
gpdd_d$bestR2=select(gpdd_d,R3a,R4a,R5a) %>% apply(1,max) #R2 for abundance
gpdd_d$bestR2m=select(gpdd_d,R3m,R4m,R5m) %>% apply(1,max) #R2 for growth rate
gpdd_results$modelresultsbest=cbind(select(gpdd_results, modelresults3, modelresults4, modelresults5),gpdd_d$bestmodel) %>% apply(1, function(x) {m=as.numeric(x["gpdd_d$bestmodel"]); x[m][[1]]})
#get E, tau, theta values from best model
gpdd_d$E=map_dbl(gpdd_results$modelresultsbest, ~.x$modelstats$E)
gpdd_d$tau=map_dbl(gpdd_results$modelresultsbest, ~.x$modelstats$tau)
gpdd_d$theta=map_dbl(gpdd_results$modelresultsbest, ~.x$modelstats$theta)
#get best model form
gpdd_d$modelform=map_chr(gpdd_results$modelresultsbest, ~.x$form)

#get jacobian matrices
gpdd_results$jacobians=map(gpdd_results$modelresultsbest, getJacobians)

#get LE estimates
gpdd_results$JLE=map2(gpdd_results$modelresultsbest, gpdd_results$jacobians, LEshift)
gpdd_d$minmean=map_dbl(gpdd_results$JLE, ~.x$minmean) #LE mean
gpdd_d$minci=map_dbl(gpdd_results$JLE, ~.x$minci) #LE lower confidence bound (*this is the LE estimate to use!*)
gpdd_d$mincisign=ifelse(gpdd_d$minci>0.01, "chaotic", "not chaotic")
length(which(gpdd_d$minci>0.01))/length(which(!is.na(gpdd_d$minci)))
#convert LE estimate to common timescale
gpdd_d$minci_mo=gpdd_d$minci/timescale_mo(gpdd_d$SamplingInterval, 1)
gpdd_d$minci_gen=gpdd_d$minci_mo*gpdd_d$MinAge_mo
gpdd_d$minmean_mo=gpdd_d$minmean/timescale_mo(gpdd_d$SamplingInterval, 1)
gpdd_d$minmean_gen=gpdd_d$minmean_mo*gpdd_d$MinAge_mo

#predictability of time series (abundance, growth rate, both, neither)
predthreshold=0.2
gpdd_d$predictable_ag=ifelse(gpdd_d$bestR2>predthreshold & gpdd_d$bestR2m>predthreshold, "ag",
                              ifelse(gpdd_d$bestR2>predthreshold & gpdd_d$bestR2m<=predthreshold, "a",
                                     ifelse(gpdd_d$bestR2<=predthreshold & gpdd_d$bestR2m>predthreshold, "g", "none")))

#fix E to 1 for best model, recompute LE ####
gpdd_results$LE1d=map2(gpdd_d$data_rescale, gpdd_d$bestmodel+2, LE1d, y="PopRescale")
gpdd_d$minci1d=map_dbl(gpdd_results$LE1d, ~.x$JLE$minci)
gpdd_d$minmean1d=map_dbl(gpdd_results$LE1d, ~.x$JLE$minmean)
gpdd_d$mincisign1d=ifelse(gpdd_d$minci1d>0.01, "chaotic", "not chaotic")
length(which(gpdd_d$minci1d>0.01))/length(which(!is.na(gpdd_d$minci1d)))

#shorten time series, recompute LE ####
gpdd1=select(gpdd_d, MainID:Notes, monotonicR2, data_rescale, bestmodel:predictable_ag) %>% 
  mutate(tslengthcat="full") %>% filter(mincisign=="chaotic")
gpdd2=gpdd1 %>% 
  mutate(datasetlength=ifelse(datasetlength==30,NA,ifelse(datasetlength/2<=30,30,ceiling(datasetlength/2))),
         tslengthcat="short") %>% filter(!is.na(datasetlength))
gpdd3=gpdd2 %>% 
  mutate(datasetlength=ifelse(datasetlength==30,NA,ifelse(datasetlength/2<=30,30,ceiling(datasetlength/2)))) %>% 
  filter(!is.na(datasetlength))
gpdd4=gpdd3 %>% 
  mutate(datasetlength=ifelse(datasetlength==30,NA,ifelse(datasetlength/2<=30,30,ceiling(datasetlength/2)))) %>% 
  filter(!is.na(datasetlength))
gpdd_short=rbind(gpdd2, gpdd3, gpdd4) %>% 
  mutate(timescale_MinAge=datasetlength/timescale_mo(SamplingInterval, MinAge_mo),
         timescale_Lifespan=datasetlength/timescale_mo(SamplingInterval, Lifespan_mo))
shortents=function(datasetlength, data) {
  data[(nrow(data)-datasetlength+1):nrow(data),]
}
gpdd_short$data_rescale=map2(gpdd_short$datasetlength, gpdd_short$data_rescale, shortents)

#run models
gpdd_short_results=select(gpdd_short, MainID, datasetlength)
gpdd_short_results$modelresults3=map(gpdd_short$data_rescale, smap_model_options, y="PopRescale", model=3)
gpdd_short_results$modelresults4=map(gpdd_short$data_rescale, smap_model_options, y="PopRescale", model=4)
gpdd_short_results$modelresults5=map(gpdd_short$data_rescale, smap_model_options, y="PopRescale", model=5)

gpdd_short_results$R3m=map_dbl(gpdd_short_results$modelresults3, ~.x$modelstats$R2model)
gpdd_short_results$R4m=map_dbl(gpdd_short_results$modelresults4, ~.x$modelstats$R2model)
gpdd_short_results$R5m=map_dbl(gpdd_short_results$modelresults5, ~.x$modelstats$R2model)
gpdd_short_results$R3a=map_dbl(gpdd_short_results$modelresults3, ~.x$modelstats$R2abund)
gpdd_short_results$R4a=map_dbl(gpdd_short_results$modelresults4, ~.x$modelstats$R2abund)
gpdd_short_results$R5a=map_dbl(gpdd_short_results$modelresults5, ~.x$modelstats$R2abund)

#best model
gpdd_short$bestmodel=select(gpdd_short_results,R3a,R4a,R5a) %>% apply(1,which.max)
gpdd_short$bestR2=select(gpdd_short_results,R3a,R4a,R5a) %>% apply(1,max)
gpdd_short$bestR2m=select(gpdd_short_results,R3m,R4m,R5m) %>% apply(1,max)
gpdd_short_results$modelresultsbest=cbind(select(gpdd_short_results, modelresults3, modelresults4, modelresults5),gpdd_short$bestmodel) %>% apply(1, function(x) {m=as.numeric(x["gpdd_short$bestmodel"]); x[m][[1]]})
gpdd_short$E=map_dbl(gpdd_short_results$modelresultsbest, ~.x$modelstats$E)
gpdd_short$tau=map_dbl(gpdd_short_results$modelresultsbest, ~.x$modelstats$tau)
gpdd_short$theta=map_dbl(gpdd_short_results$modelresultsbest, ~.x$modelstats$theta)
gpdd_short$modelform=map_chr(gpdd_short_results$modelresultsbest, ~.x$form)

#jacobians and LEs
gpdd_short_results$jacobians=map(gpdd_short_results$modelresultsbest, getJacobians)
gpdd_short_results$JLE=map2(gpdd_short_results$modelresultsbest, gpdd_short_results$jacobians, LEshift)
gpdd_short$minci=map_dbl(gpdd_short_results$JLE, ~.x$minci)
gpdd_short$minmean=map_dbl(gpdd_short_results$JLE, ~.x$minmean)
gpdd_short$mincisign=ifelse(gpdd_short$minci>0.01, "chaotic", "not chaotic")
#convert to common timescale
gpdd_short$minci_mo=gpdd_short$minci/timescale_mo(gpdd_short$SamplingInterval, 1)
gpdd_short$minci_gen=gpdd_short$minci_mo*gpdd_short$MinAge_mo
gpdd_short$minmean_mo=gpdd_short$minmean/timescale_mo(gpdd_short$SamplingInterval, 1)
gpdd_short$minmean_gen=gpdd_short$minmean_mo*gpdd_short$MinAge_mo

#predictability
gpdd_short$predictable_ag=ifelse(gpdd_short$bestR2>predthreshold & gpdd_short$bestR2m>predthreshold, "ag",
                             ifelse(gpdd_short$bestR2>predthreshold & gpdd_short$bestR2m<=predthreshold, "a",
                                    ifelse(gpdd_short$bestR2<=predthreshold & gpdd_short$bestR2m>predthreshold, "g", "none")))

gpdd_combo=rbind(gpdd1,gpdd_short)

# sensitivity to changes in E and tau ####
gpdd_d$modnum=case_when(gpdd_d$modelform =="ut-ut" ~ 1,
                        gpdd_d$modelform =="log-log" ~ 2,
                        gpdd_d$modelform =="fd-ut" ~ 3,
                        gpdd_d$modelform =="gr-ut" ~ 4,
                        gpdd_d$modelform =="gr-log" ~ 5)
gpdd_d$Eplus=ifelse(gpdd_d$E+1>6,6,gpdd_d$E+1)
gpdd_d$Eminus=ifelse(gpdd_d$E-1<1,1,gpdd_d$E-1)
gpdd_d$tauplus=ifelse(gpdd_d$tau+1>6,6,gpdd_d$tau+1)
gpdd_d$tauminus=ifelse(gpdd_d$tau-1<1,1,gpdd_d$tau-1)

gpdd_results$sens1=pmap(list(data=gpdd_d$data_rescale, model=gpdd_d$modnum,
                             Efix=gpdd_d$Eplus,taufix=gpdd_d$tau),LEfix, y="PopRescale")
gpdd_results$sens2=pmap(list(data=gpdd_d$data_rescale, model=gpdd_d$modnum,
                             Efix=gpdd_d$Eminus,taufix=gpdd_d$tau),LEfix, y="PopRescale")
gpdd_results$sens3=pmap(list(data=gpdd_d$data_rescale, model=gpdd_d$modnum,
                             Efix=gpdd_d$E,taufix=gpdd_d$tauplus),LEfix, y="PopRescale")
gpdd_results$sens4=pmap(list(data=gpdd_d$data_rescale, model=gpdd_d$modnum,
                             Efix=gpdd_d$E,taufix=gpdd_d$tauminus),LEfix, y="PopRescale")

gpdd_d$sens1LEmin=map_dbl(gpdd_results$sens1, ~.x$JLE$minci)
gpdd_d$sens2LEmin=map_dbl(gpdd_results$sens2, ~.x$JLE$minci)
gpdd_d$sens3LEmin=map_dbl(gpdd_results$sens3, ~.x$JLE$minci)
gpdd_d$sens4LEmin=map_dbl(gpdd_results$sens4, ~.x$JLE$minci)

gpdd_d$sens1LEsign=ifelse(gpdd_d$sens1LEmin>0.01, "chaotic", "not chaotic")
gpdd_d$sens2LEsign=ifelse(gpdd_d$sens2LEmin>0.01, "chaotic", "not chaotic")
gpdd_d$sens3LEsign=ifelse(gpdd_d$sens3LEmin>0.01, "chaotic", "not chaotic")
gpdd_d$sens4LEsign=ifelse(gpdd_d$sens4LEmin>0.01, "chaotic", "not chaotic")

gpdd_d$sens1R2=map_dbl(gpdd_results$sens1, ~.x$modelresults$modelstats$R2abund)
gpdd_d$sens2R2=map_dbl(gpdd_results$sens2, ~.x$modelresults$modelstats$R2abund)
gpdd_d$sens3R2=map_dbl(gpdd_results$sens3, ~.x$modelresults$modelstats$R2abund)
gpdd_d$sens4R2=map_dbl(gpdd_results$sens4, ~.x$modelresults$modelstats$R2abund)

sensresults=data.frame(model=c("best","E+1","E-1","tau+1","tau-1"), propchaotic=NA, meanR2=NA, medianR2=NA)

sensresults$propchaotic[1]=length(which(gpdd_d$mincisign=="chaotic"))/length(which(!is.na(gpdd_d$mincisign)))
sensresults$propchaotic[2]=length(which(gpdd_d$sens1LEsign=="chaotic"))/length(which(!is.na(gpdd_d$mincisign)))
sensresults$propchaotic[3]=length(which(gpdd_d$sens2LEsign=="chaotic"))/length(which(!is.na(gpdd_d$mincisign)))
sensresults$propchaotic[4]=length(which(gpdd_d$sens3LEsign=="chaotic"))/length(which(!is.na(gpdd_d$mincisign)))
sensresults$propchaotic[5]=length(which(gpdd_d$sens4LEsign=="chaotic"))/length(which(!is.na(gpdd_d$mincisign)))

sensresults$meanR2[1]=mean(gpdd_d$bestR2)
sensresults$meanR2[2]=mean(gpdd_d$sens1R2)
sensresults$meanR2[3]=mean(gpdd_d$sens2R2)
sensresults$meanR2[4]=mean(gpdd_d$sens3R2)
sensresults$meanR2[5]=mean(gpdd_d$sens4R2)

sensresults$medianR2[1]=median(gpdd_d$bestR2)
sensresults$medianR2[2]=median(gpdd_d$sens1R2)
sensresults$medianR2[3]=median(gpdd_d$sens2R2)
sensresults$medianR2[4]=median(gpdd_d$sens3R2)
sensresults$medianR2[5]=median(gpdd_d$sens4R2)

sensresults[,2:4]=round(sensresults[,2:4],2)

#### Export Results ####

#save results
save(gpdd_d, gpdd_results, gpdd_combo, gpdd_short, gpdd_short_results, file = "./data/gpdd_results_update3.Rdata")

#export E and tau for use in other analyses
exportEtau=select(gpdd_d, MainID, E, tau)
write.csv(exportEtau, "./data/gpdd_Etau_smap.csv", row.names = F)

#export main results
#change names of some variables, so more intuitive
exportres=select(gpdd_d, MainID, R2abund=bestR2, R2gr=bestR2m, predictable_ag, modelform, E, tau, theta, 
                 LEmean=minmean, LEmean_mo=minmean_mo, LEmean_gen=minmean_gen, 
                 LEmin=minci, LEmin_mo=minci_mo, LEmin_gen=minci_gen, LEclass=mincisign, 
                 LEmin1d=minci1d, LEmean1d=minmean1d, LEclass1d=mincisign1d)
write.csv(exportres, "./data/gpdd_results_smap.csv", row.names = F)
#export results with shortened time series
exportres2=select(gpdd_combo, MainID, datasetlength, tslengthcat, timescale_MinAge, MinAge_mo, Mass_g, R2abund=bestR2, R2gr=bestR2m, predictable_ag, E, tau, theta, 
                 LEmean=minmean, LEmin=minci, LEmin_mo=minci_mo, LEmin_gen=minci_gen, LEclass=mincisign)
write.csv(exportres2, "./data/gpdd_results_truncation_smap.csv", row.names = F)
