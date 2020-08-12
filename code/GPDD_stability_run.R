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

#fit models
gpdd_results=select(gpdd_d, MainID)
#ut-ut
gpdd_results$modelresults1=map(gpdd_d$data_rescale, smap_model, y="PopRescale", ylog=F)
#log-log
gpdd_results$modelresults2=map(gpdd_d$data_rescale, smap_model, y="PopRescale_log", ylog=T)
#fd-ut
gpdd_results$modelresults3=map(gpdd_d$data_rescale, smap_model_gr_fd, fd="PopRescale_fd", y="PopRescale", ylog=F)
#gr-ut
gpdd_results$modelresults4=map(gpdd_d$data_rescale, smap_model_gr_fd, gr="PopRescale_gr", y="PopRescale", ylog=F)
#gr-log
gpdd_results$modelresults5=map(gpdd_d$data_rescale, smap_model_gr_fd, gr="PopRescale_gr", y="PopRescale_log", ylog=T)

#pull out R2 values for each model
gpdd_d$R1m=map_dbl(gpdd_results$modelresults1, ~.x$modelstats$R2model)
gpdd_d$R1a=map_dbl(gpdd_results$modelresults1, ~.x$modelstats$R2abund)
gpdd_d$R2m=map_dbl(gpdd_results$modelresults2, ~.x$modelstats$R2model)
gpdd_d$R2a=map_dbl(gpdd_results$modelresults2, ~.x$modelstats$R2abund)
gpdd_d$R3m=map_dbl(gpdd_results$modelresults3, ~.x$modelstats$R2model)
gpdd_d$R3a=map_dbl(gpdd_results$modelresults3, ~.x$modelstats$R2abund)
gpdd_d$R4m=map_dbl(gpdd_results$modelresults4, ~.x$modelstats$R2model)
gpdd_d$R4a=map_dbl(gpdd_results$modelresults4, ~.x$modelstats$R2abund)
gpdd_d$R5m=map_dbl(gpdd_results$modelresults5, ~.x$modelstats$R2model)
gpdd_d$R5a=map_dbl(gpdd_results$modelresults5, ~.x$modelstats$R2abund)
# gpdd_d=cbind(gpdd_d,map_df(gpdd_d$modelresults1, ~.x$modelstats))

#best model
gpdd_d$bestmodel=select(gpdd_d,R1a,R2a,R3a,R4a,R5a) %>% apply(1,which.max)
gpdd_results$modelresultsbest=cbind(select(gpdd_results, modelresults1:modelresults5),gpdd_d$bestmodel) %>% apply(1, function(x) {m=as.numeric(x["gpdd_d$bestmodel"]); x[m][[1]]})

#predictability
predthreshold=0.2
gpdd_d$predictable=select(gpdd_d,R1a,R2a,R3a,R4a,R5a) %>% apply(1,FUN = function(x) ifelse(any(x>predthreshold), "yes","no"))
gpdd_d$predictable_log=ifelse(gpdd_d$R2m>predthreshold, "yes","no")
gpdd_d$predictable_gr=select(gpdd_d,R3m,R4m,R5m) %>% apply(1,FUN = function(x) ifelse(any(x>predthreshold), "yes","no"))
gpdd_d$predictable_ag=ifelse(gpdd_d$predictable=="yes" & gpdd_d$predictable_gr=="yes", "ag",
                              ifelse(gpdd_d$predictable=="yes" & gpdd_d$predictable_gr=="no", "a",
                                     ifelse(gpdd_d$predictable=="no" & gpdd_d$predictable_gr=="yes", "g", "none")))

#get jacobians
gpdd_results$jacobians1=map(gpdd_results$modelresults1, getJacobians)
gpdd_results$jacobians2=map(gpdd_results$modelresults2, getJacobians)
gpdd_results$jacobians3=map(gpdd_results$modelresults3, getJacobians)
gpdd_results$jacobians4=map(gpdd_results$modelresults4, getJacobians)
gpdd_results$jacobians5=map(gpdd_results$modelresults5, getJacobians)
#get stability
gpdd_results$stability1=map(gpdd_results$jacobians1, getStability)
gpdd_results$stability2=map(gpdd_results$jacobians2, getStability)
gpdd_results$stability3=map(gpdd_results$jacobians3, getStability)
gpdd_results$stability4=map(gpdd_results$jacobians4, getStability)
gpdd_results$stability5=map(gpdd_results$jacobians5, getStability)

gpdd_d$gle1=map_dbl(gpdd_results$stability1, ~.x$gle)
gpdd_d$gle2=map_dbl(gpdd_results$stability2, ~.x$gle)
gpdd_d$gle3=map_dbl(gpdd_results$stability3, ~.x$gle)
gpdd_d$gle4=map_dbl(gpdd_results$stability4, ~.x$gle)
gpdd_d$gle5=map_dbl(gpdd_results$stability5, ~.x$gle)
gpdd_d$lle_avg1=map_dbl(gpdd_results$stability1, ~.x$lle_avg)
gpdd_d$lle_avg2=map_dbl(gpdd_results$stability2, ~.x$lle_avg)
gpdd_d$lle_avg3=map_dbl(gpdd_results$stability3, ~.x$lle_avg)
gpdd_d$lle_avg4=map_dbl(gpdd_results$stability4, ~.x$lle_avg)
gpdd_d$lle_avg5=map_dbl(gpdd_results$stability5, ~.x$lle_avg)
gpdd_d$lle_pp1=map_dbl(gpdd_results$stability1, ~.x$lle_pp)
gpdd_d$lle_pp2=map_dbl(gpdd_results$stability2, ~.x$lle_pp)
gpdd_d$lle_pp3=map_dbl(gpdd_results$stability3, ~.x$lle_pp)
gpdd_d$lle_pp4=map_dbl(gpdd_results$stability4, ~.x$lle_pp)
gpdd_d$lle_pp5=map_dbl(gpdd_results$stability5, ~.x$lle_pp)

#best stability estimates
gpdd_d$glebest=apply(gpdd_d, 1, FUN=function(x) {unlist(x[paste0("gle",x$bestmodel)])})
gpdd_d$lle_avgbest= apply(gpdd_d, 1, FUN=function(x) {unlist(x[paste0("lle_avg",x$bestmodel)])})
gpdd_d$lle_ppbest= apply(gpdd_d, 1, FUN=function(x) {unlist(x[paste0("lle_pp",x$bestmodel)])})
gpdd_d$glesign=ifelse(gpdd_d$glebest>0, "positive", "negative")
gpdd_d$lle_avgsign=ifelse(gpdd_d$lle_avgbest>0, "positive", "negative")
gpdd_d$Ebest=map_dbl(gpdd_results$modelresultsbest, ~.x$modelstats$E)
gpdd_d$thetabest=map_dbl(gpdd_results$modelresultsbest, ~.x$modelstats$theta)

#convert to common timescale
gpdd_d$gle_mo=gpdd_d$glebest/timescale_mo(gpdd_d$SamplingInterval, 1)
gpdd_d$gle_gen=gpdd_d$gle_mo*gpdd_d$MinAge_mo

#regression method
gpdd_results$regLE=map(gpdd_d$data_rescale, regLE, y="PopRescale")
gpdd_d$LEreg=map_dbl(gpdd_results$regLE, ~.x$LEreg)
gpdd_d$LEreg_se=map_dbl(gpdd_results$regLE, ~.x$LEreg_se)

#fix E to 1 for best model
gpdd_results$LE1d=map2(gpdd_d$data_rescale, gpdd_d$bestmodel, LE1d)
gpdd_d$glebest1d=map_dbl(gpdd_results$LE1d, ~.x$gle)
gpdd_d$glesign1d=ifelse(gpdd_d$glebest1d>0, "positive", "negative")
gpdd_d$gle1d_mo=gpdd_d$glebest1d/timescale_mo(gpdd_d$SamplingInterval, 1)
gpdd_d$gle1d_gen=gpdd_d$gle1d_mo*gpdd_d$MinAge_mo

#prediction slope for best model



save(gpdd_d, gpdd_results, file = "./Data/gpdd_results.Rdata")

# #increase threshold
# predthreshold=0.2
# gpdd_d$predictable=select(gpdd_d,R1a,R2a,R3a,R4a,R5a) %>% apply(1,FUN = function(x) ifelse(any(x>predthreshold), "yes","no"))
# gpdd_d$predictable_gr=select(gpdd_d,R3m,R4m,R5m) %>% apply(1,FUN = function(x) ifelse(any(x>predthreshold), "yes","no"))
# gpdd_d$predictable_ag=ifelse(gpdd_d$predictable=="yes" & gpdd_d2$predictable_gr=="yes", "ag",
#                              ifelse(gpdd_d2$predictable=="yes" & gpdd_d2$predictable_gr=="no", "a",
#                                     ifelse(gpdd_d2$predictable=="no" & gpdd_d2$predictable_gr=="yes", "g", "none")))
# gpdd_d3=filter(gpdd_d, predictable_ag!="none")
# gpdd_results3=filter(gpdd_results, MainID %in% unique(gpdd_d3$MainID))
