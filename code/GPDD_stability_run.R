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

gpdd_results=select(gpdd_d, MainID)

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

#get best model
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
gpdd_results$LEshift=pmap(list(gpdd_results$modelresultsbest, gpdd_results$jacobians, gpdd_d$SamplingInterval), LEshift)
#LEs of best model
gpdd_d$gle=map_dbl(gpdd_results$stability, ~.x$gle)
gpdd_d$minci=map_dbl(gpdd_results$LEshift, ~.x$minci)
gpdd_d$minmean=map_dbl(gpdd_results$LEshift, ~.x$minmean)
gpdd_d$lle_pp=map_dbl(gpdd_results$stability, ~.x$lle_pp)
#LE signs
gpdd_d$glesign=ifelse(gpdd_d$gle>0.01, "chaotic", "not chaotic")
gpdd_d$mincisign=ifelse(gpdd_d$minci>0.01, "chaotic", "not chaotic")
gpdd_d$minmeansign=ifelse(gpdd_d$minmean>0.01, "chaotic", "not chaotic")

length(which(gpdd_d$gle>0.01))/length(which(!is.na(gpdd_d$gle)))
length(which(gpdd_d$minci>0.01))/length(which(!is.na(gpdd_d$minci)))
length(which(gpdd_d$minmean>0.01))/length(which(!is.na(gpdd_d$minmean)))

#predictability
predthreshold=0.2
gpdd_d$predictable=ifelse(gpdd_d$bestR2>predthreshold, "yes","no")
gpdd_d$predictable_gr=ifelse(gpdd_d$bestR2m>predthreshold, "yes","no")
gpdd_d$predictable_ag=ifelse(gpdd_d$predictable=="yes" & gpdd_d$predictable_gr=="yes", "ag",
                              ifelse(gpdd_d$predictable=="yes" & gpdd_d$predictable_gr=="no", "a",
                                     ifelse(gpdd_d$predictable=="no" & gpdd_d$predictable_gr=="yes", "g", "none")))

# #verify 1&3 and 2&5 are the same
# gpdd_results2$jacobians1=map(gpdd_results$modelresults1, getJacobians)
# gpdd_results2$jacobians2=map(gpdd_results$modelresults2, getJacobians)
# gpdd_results2$jacobians3=map(gpdd_results$modelresults3, getJacobians)
# gpdd_results2$jacobians4=map(gpdd_results$modelresults4, getJacobians)
# gpdd_results2$jacobians5=map(gpdd_results$modelresults5, getJacobians)
# gpdd_results2$stability1=map2(gpdd_results$modelresults1, gpdd_results2$jacobians1, getStability)
# gpdd_results2$stability2=map2(gpdd_results$modelresults2, gpdd_results2$jacobians2, getStability)
# gpdd_results2$stability3=map2(gpdd_results$modelresults3, gpdd_results2$jacobians3, getStability)
# gpdd_results2$stability4=map2(gpdd_results$modelresults4, gpdd_results2$jacobians4, getStability)
# gpdd_results2$stability5=map2(gpdd_results$modelresults5, gpdd_results2$jacobians5, getStability)
# gpdd_d2$gle1=map_dbl(gpdd_results2$stability1, ~.x$gle)
# gpdd_d2$gle2=map_dbl(gpdd_results2$stability2, ~.x$gle)
# gpdd_d2$gle3=map_dbl(gpdd_results2$stability3, ~.x$gle)
# gpdd_d2$gle4=map_dbl(gpdd_results2$stability4, ~.x$gle)
# gpdd_d2$gle5=map_dbl(gpdd_results2$stability5, ~.x$gle)

#convert to common timescale
gpdd_d$gle_mo=gpdd_d$gle/timescale_mo(gpdd_d$SamplingInterval, 1)
gpdd_d$gle_gen=gpdd_d$gle_mo*gpdd_d$MinAge_mo
gpdd_d$minci_mo=gpdd_d$minci/timescale_mo(gpdd_d$SamplingInterval, 1)
gpdd_d$minci_gen=gpdd_d$minci_mo*gpdd_d$MinAge_mo
gpdd_d$minmean_mo=gpdd_d$minci/timescale_mo(gpdd_d$SamplingInterval, 1)
gpdd_d$minmean_gen=gpdd_d$minmean_mo*gpdd_d$MinAge_mo

#fix E to 1 for best model
#tau unconstrained
gpdd_results$LE1d=map2(gpdd_d$data_rescale, gpdd_d$bestmodel+2, LE1d, y="PopRescale")
gpdd_results$jacobians1d=map(gpdd_results$LE1d, ~.x$jacobians)
gpdd_results$modelresults1d=map(gpdd_results$LE1d, ~.x$modelresults)
gpdd_results$LEshift1d=pmap(list(gpdd_results$modelresults1d, gpdd_results$jacobians1d, gpdd_d$SamplingInterval), LEshift)

gpdd_d$gle1d=map_dbl(gpdd_results$LE1d, ~.x$stability$gle)
gpdd_d$minci1d=map_dbl(gpdd_results$LEshift1d, ~.x$minci)
gpdd_d$minmean1d=map_dbl(gpdd_results$LEshift1d, ~.x$minmean)
gpdd_d$mincisign1d=ifelse(gpdd_d$minci1d>0.01, "chaotic", "not chaotic")

length(which(gpdd_d$gle1d>0.01))/length(which(!is.na(gpdd_d$gle1d)))
length(which(gpdd_d$minci1d>0.01))/length(which(!is.na(gpdd_d$minci1d)))

#use shortened time series ####

gpdd1=select(gpdd_d, MainID:Notes, monotonicR2, data_rescale, bestmodel:minmean_gen) %>% 
  mutate(tslengthcat="full")
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
  mutate(timescale_MinAge=timescale_ratio(SamplingInterval, datasetlength, MinAge_mo),
         timescale_Lifespan=timescale_ratio(SamplingInterval, datasetlength, Lifespan_mo))
shortents=function(datasetlength, data) {
  data[(nrow(data)-datasetlength+1):nrow(data),]
}
gpdd_short$data_rescale=map2(gpdd_short$datasetlength, gpdd_short$data_rescale, shortents)

#run models
gpdd_short_results=select(gpdd_short, MainID, datasetlength)
#fd-ut
gpdd_short_results$modelresults3=map(gpdd_short$data_rescale, smap_model_options, y="PopRescale", model=3)
#gr-ut
gpdd_short_results$modelresults4=map(gpdd_short$data_rescale, smap_model_options, y="PopRescale", model=4)
#gr-log
gpdd_short_results$modelresults5=map(gpdd_short$data_rescale, smap_model_options, y="PopRescale", model=5)

gpdd_short_results$R3m=map_dbl(gpdd_short_results$modelresults3, ~.x$modelstats$R2model)
gpdd_short_results$R4m=map_dbl(gpdd_short_results$modelresults4, ~.x$modelstats$R2model)
gpdd_short_results$R5m=map_dbl(gpdd_short_results$modelresults5, ~.x$modelstats$R2model)
gpdd_short_results$R3a=map_dbl(gpdd_short_results$modelresults3, ~.x$modelstats$R2abund)
gpdd_short_results$R4a=map_dbl(gpdd_short_results$modelresults4, ~.x$modelstats$R2abund)
gpdd_short_results$R5a=map_dbl(gpdd_short_results$modelresults5, ~.x$modelstats$R2abund)
# gpdd_short=cbind(gpdd_short,map_df(gpdd_short$modelresults1, ~.x$modelstats))

#best model
gpdd_short$bestmodel=select(gpdd_short_results,R3a,R4a,R5a) %>% apply(1,which.max)
gpdd_short$bestR2=select(gpdd_short_results,R3a,R4a,R5a) %>% apply(1,max)
gpdd_short$bestR2m=select(gpdd_short_results,R3m,R4m,R5m) %>% apply(1,max)
gpdd_short_results$modelresultsbest=cbind(select(gpdd_short_results, modelresults3, modelresults4, modelresults5),gpdd_short$bestmodel) %>% apply(1, function(x) {m=as.numeric(x["gpdd_short$bestmodel"]); x[m][[1]]})
gpdd_short$E=map_dbl(gpdd_short_results$modelresultsbest, ~.x$modelstats$E)
gpdd_short$tau=map_dbl(gpdd_short_results$modelresultsbest, ~.x$modelstats$tau)
gpdd_short$theta=map_dbl(gpdd_short_results$modelresultsbest, ~.x$modelstats$theta)

#jacobians and stability of best model
gpdd_short_results$jacobians=map(gpdd_short_results$modelresultsbest, getJacobians)
gpdd_short_results$stability=map2(gpdd_short_results$modelresultsbest, gpdd_short_results$jacobians, getStability)
gpdd_short_results$LEshift=pmap(list(gpdd_short_results$modelresultsbest, gpdd_short_results$jacobians, gpdd_short$SamplingInterval), LEshift)

#LEs of best model
gpdd_short$gle=map_dbl(gpdd_short_results$stability, ~.x$gle)
gpdd_short$minci=map_dbl(gpdd_short_results$LEshift, ~.x$minci)
gpdd_short$minmean=map_dbl(gpdd_short_results$LEshift, ~.x$minmean)
gpdd_short$lle_pp=map_dbl(gpdd_short_results$stability, ~.x$lle_pp)
#LE signs
gpdd_short$glesign=ifelse(gpdd_short$gle>0.01, "chaotic", "not chaotic")
gpdd_short$mincisign=ifelse(gpdd_short$minci>0.01, "chaotic", "not chaotic")
gpdd_short$minmeansign=ifelse(gpdd_short$minmean>0.01, "chaotic", "not chaotic")
#predictability
predthreshold=0.2
gpdd_short$predictable=ifelse(gpdd_short$bestR2>predthreshold, "yes","no")
gpdd_short$predictable_gr=ifelse(gpdd_short$bestR2m>predthreshold, "yes","no")
gpdd_short$predictable_ag=ifelse(gpdd_short$predictable=="yes" & gpdd_short$predictable_gr=="yes", "ag",
                             ifelse(gpdd_short$predictable=="yes" & gpdd_short$predictable_gr=="no", "a",
                                    ifelse(gpdd_short$predictable=="no" & gpdd_short$predictable_gr=="yes", "g", "none")))
#convert to common timescale
gpdd_short$gle_mo=gpdd_short$gle/timescale_mo(gpdd_short$SamplingInterval, 1)
gpdd_short$gle_gen=gpdd_short$gle_mo*gpdd_short$MinAge_mo
gpdd_short$minci_mo=gpdd_short$minci/timescale_mo(gpdd_short$SamplingInterval, 1)
gpdd_short$minci_gen=gpdd_short$minci_mo*gpdd_short$MinAge_mo
gpdd_short$minmean_mo=gpdd_short$minci/timescale_mo(gpdd_short$SamplingInterval, 1)
gpdd_short$minmean_gen=gpdd_short$minmean_mo*gpdd_short$MinAge_mo

gpdd_combo=rbind(gpdd1,gpdd_short)

#regression method
gpdd_results$regLE=map2(gpdd_d$data_rescale, gpdd_results$modelresultsbest, regLE, y="PopRescale")
gpdd_d$LEreg=map_dbl(gpdd_results$regLE, ~.x$LEreg)
gpdd_d$LEreg_se=map_dbl(gpdd_results$regLE, ~.x$LEreg_se)
gpdd_d$LEregmin=gpdd_d$LEreg-1.96*gpdd_d$LEreg_se

#sibly method
gpdd_d$LEsibly=map_dbl(gpdd_d$data_rescale, SiblymodelLE, y="PopRescale")

#get variance on LE
gpdd_d$varmin=map_dbl(gpdd_results$LEshift, ~.x$varmin)
gpdd_d$varmin_mo=map_dbl(gpdd_results$LEshift, ~.x$varmin_mo)

#get variance intercept
gpdd_results$LEsaturation=pmap(list(gpdd_results$modelresultsbest, gpdd_results$jacobians, gpdd_d$SamplingInterval), LEsaturation)
gpdd_d$logerror=map_dbl(gpdd_results$LEsaturation, ~.x$logerror)
gpdd_d$slope=map_dbl(gpdd_results$LEsaturation, ~.x$slope)

#get best model form
gpdd_d$modelform=map_chr(gpdd_results$modelresultsbest, ~.x$form)

#repeat with sqrt transform
gpdd_results$modelresultssqrt=map(gpdd_d$data_sqrt, smap_model_options, y="PopSqrt", model=3)
gpdd_results$jacobianssqrt=map(gpdd_results$modelresultssqrt, getJacobians)
gpdd_results$LEshiftsqrt=pmap(list(gpdd_results$modelresultssqrt, gpdd_results$jacobianssqrt, gpdd_d$SamplingInterval), LEshift)
#LEs of best model
gpdd_d$LEmin_sqrt=map_dbl(gpdd_results$LEshiftsqrt, ~.x$minci)
#LE signs
gpdd_d$LEclass_sqrt=ifelse(gpdd_d$LEmin_sqrt>0.01, "chaotic", "not chaotic")
gpdd_d$R2sqrt=map_dbl(gpdd_results$modelresultssqrt, ~.x$modelstats$R2model)
#
gpdd_results$jacobians3=map(gpdd_results$modelresults3, getJacobians)
gpdd_results$LEshift3=pmap(list(gpdd_results$modelresults3, gpdd_results$jacobians3, gpdd_d$SamplingInterval), LEshift)
gpdd_d$LEmin_3=map_dbl(gpdd_results$LEshift3, ~.x$minci)
#
gpdd_results$modelresultsut=map(gpdd_d$data_sqrt, smap_model_options, y="Pop", model=3)
gpdd_results$jacobiansut=map(gpdd_results$modelresultsut, getJacobians)
gpdd_results$LEshiftut=pmap(list(gpdd_results$modelresultsut, gpdd_results$jacobiansut, gpdd_d$SamplingInterval), LEshift)
gpdd_d$LEmin_ut=map_dbl(gpdd_results$LEshiftut, ~.x$minci)

#growth rates
gpdd_d$meangr=map_dbl(gpdd_d$data_rescale, ~mean(.x$PopRescale_gr, na.rm=T))
gpdd_d$maxgr=map_dbl(gpdd_d$data_rescale, ~max(.x$PopRescale_gr, na.rm=T))
gpdd_d$maxgr_mo=gpdd_d$maxgr/timescale_mo(gpdd_d$SamplingInterval,1)

#save
save(gpdd_d, gpdd_results, gpdd_combo, gpdd_short, gpdd_short_results, file = "./data/gpdd_results_update2.Rdata")

#export E and tau for other analyses
exportEtau=select(gpdd_d, MainID, E, tau)
write.csv(exportEtau, "./data/gpdd_Etau_smap.csv", row.names = F)

#export results
exportres=select(gpdd_d, MainID, R2abund=bestR2, R2gr=bestR2m, predictable_ag, modelform, E, tau, theta, 
                 LEmean=minmean, LEmin=minci, LEmin_mo=minci_mo, LEmin_gen=minci_gen, LEclass=mincisign, 
                 LEmean1d=minmean1d, LEmin1d=minci1d, LEclass1d=mincisign1d, LLE_proppos=lle_pp, LEreg, LEreg_se, LEregmin, 
                 LEsibly, LEvar=varmin, LEvar_mo=varmin_mo, logerror, slope, maxgr_mo)
write.csv(exportres, "./data/gpdd_results_smap.csv", row.names = F)
exportres2=select(gpdd_combo, MainID, datasetlength, tslengthcat, timescale_MinAge, MinAge_mo, Mass_g, R2abund=bestR2, R2gr=bestR2m, predictable_ag, E, tau, theta, 
                 LEmean=minmean, LEmin=minci, LEmin_mo=minci_mo, LEmin_gen=minci_gen, LEclass=mincisign)
write.csv(exportres2, "./data/gpdd_results_truncation_smap.csv", row.names = F)


#troubleshooting

# LEshiftres=NULL
# for(i in 1:nrow(gpdd_short_results)) {
#   LEshiftres[[i]]=LEshift(gpdd_short_results$modelresultsbest[[i]], gpdd_short_results$jacobians[[i]])
# }
# 
# #plot a time series
# plotMainID=function(ID) {
#   testplot=filter(gpdd_d, MainID==ID)
#   plot(testplot$data_rescale[[1]]$SeriesStep, testplot$data_rescale[[1]]$PopRescale, 
#        ylab="PopRescale", xlab="SeriesStep", main=paste(ID, testplot$CommonName))
#   lines(testplot$data_rescale[[1]]$SeriesStep, testplot$data_rescale[[1]]$PopRescale)
# }
# plotMainID(9953)
# plotMainID(56)

# #increase threshold
# predthreshold=0.2
# gpdd_d$predictable=select(gpdd_d,R1a,R2a,R3a,R4a,R5a) %>% apply(1,FUN = function(x) ifelse(any(x>predthreshold), "yes","no"))
# gpdd_d$predictable_gr=select(gpdd_d,R3m,R4m,R5m) %>% apply(1,FUN = function(x) ifelse(any(x>predthreshold), "yes","no"))
# gpdd_d$predictable_ag=ifelse(gpdd_d$predictable=="yes" & gpdd_d2$predictable_gr=="yes", "ag",
#                              ifelse(gpdd_d2$predictable=="yes" & gpdd_d2$predictable_gr=="no", "a",
#                                     ifelse(gpdd_d2$predictable=="no" & gpdd_d2$predictable_gr=="yes", "g", "none")))
# gpdd_d3=filter(gpdd_d, predictable_ag!="none")
# gpdd_results3=filter(gpdd_results, MainID %in% unique(gpdd_d3$MainID))
