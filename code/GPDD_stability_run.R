# GPDD stability analysis - run

library("rgpdd")
library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")
library("purrr")

source("./Code/GPDD_stability_functions.R")
load(file = "./Data/gpdd_d_v2.Rdata")
source("~/GRAD SCHOOL/R reference/ggplot themes rogers.R")

# #test dataset
# randomseries=sample(gpdd_d$MainID,30)
# #gpdd_d_TEST=gpdd_d[1:10,]
# gpdd_d_TEST=filter(gpdd_d, MainID %in% randomseries)

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

#predictability
predthreshold=0.2
gpdd_d$bestmodel=select(gpdd_d,R1a,R2a,R3a,R4a,R5a) %>% apply(1,which.max)
gpdd_d$predictable=select(gpdd_d,R1a,R2a,R3a,R4a,R5a) %>% apply(1,FUN = function(x) ifelse(any(x>predthreshold), "yes","no"))
gpdd_d$predictable_gr=select(gpdd_d,R3m,R4m,R5m) %>% apply(1,FUN = function(x) ifelse(any(x>predthreshold), "yes","no"))
gpdd_d$predictable_ag=ifelse(gpdd_d$predictable=="yes" & gpdd_d$predictable_gr=="yes", "ag",
                              ifelse(gpdd_d$predictable=="yes" & gpdd_d$predictable_gr=="no", "a",
                                     ifelse(gpdd_d$predictable=="no" & gpdd_d$predictable_gr=="yes", "g", "none")))

# gpdd_d2$predictability=ifelse(gpdd_d2$R2a>0.2 & gpdd_d2$R3m>0.2, "a g",
#                               ifelse(gpdd_d2$R2a>0.2 & gpdd_d2$R3m<=0.2, "a",
#                                      ifelse(gpdd_d2$R2a<=0.2 & gpdd_d2$R3m>0.2, "g", NA)))

#filter out unpredictable time series

#gpdd_d2=filter(gpdd_d, R2a>0.2 | R3m>0.2)
gpdd_d2=filter(gpdd_d, predictable_ag!="none")
gpdd_results2=filter(gpdd_results, MainID %in% unique(gpdd_d2$MainID))

#get jacobians
gpdd_results2$jacobians1=map(gpdd_results2$modelresults1, getJacobians)
gpdd_results2$jacobians2=map(gpdd_results2$modelresults2, getJacobians)
gpdd_results2$jacobians3=map(gpdd_results2$modelresults3, getJacobians)
gpdd_results2$jacobians4=map(gpdd_results2$modelresults4, getJacobians)
gpdd_results2$jacobians5=map(gpdd_results2$modelresults5, getJacobians)
#get stability
gpdd_results2$stability1=map(gpdd_results2$jacobians1, getStability)
gpdd_results2$stability2=map(gpdd_results2$jacobians2, getStability)
gpdd_results2$stability3=map(gpdd_results2$jacobians3, getStability)
gpdd_results2$stability4=map(gpdd_results2$jacobians4, getStability)
gpdd_results2$stability5=map(gpdd_results2$jacobians5, getStability)

gpdd_d2$lle_pp1=map_dbl(gpdd_results2$stability1, ~.x$lle_pp)
gpdd_d2$lle_pp2=map_dbl(gpdd_results2$stability2, ~.x$lle_pp)
gpdd_d2$lle_pp3=map_dbl(gpdd_results2$stability3, ~.x$lle_pp)
gpdd_d2$lle_pp4=map_dbl(gpdd_results2$stability4, ~.x$lle_pp)
gpdd_d2$lle_pp5=map_dbl(gpdd_results2$stability5, ~.x$lle_pp)
gpdd_d2$lle_avg1=map_dbl(gpdd_results2$stability1, ~.x$lle_avg)
gpdd_d2$lle_avg2=map_dbl(gpdd_results2$stability2, ~.x$lle_avg)
gpdd_d2$lle_avg3=map_dbl(gpdd_results2$stability3, ~.x$lle_avg)
gpdd_d2$lle_avg4=map_dbl(gpdd_results2$stability4, ~.x$lle_avg)
gpdd_d2$lle_avg5=map_dbl(gpdd_results2$stability5, ~.x$lle_avg)
gpdd_d2$gle1=map_dbl(gpdd_results2$stability1, ~.x$gle)
gpdd_d2$gle2=map_dbl(gpdd_results2$stability2, ~.x$gle)
gpdd_d2$gle3=map_dbl(gpdd_results2$stability3, ~.x$gle)
gpdd_d2$gle4=map_dbl(gpdd_results2$stability4, ~.x$gle)
gpdd_d2$gle5=map_dbl(gpdd_results2$stability5, ~.x$gle)

save(gpdd_d2, gpdd_results, gpdd_results2, file = "./Data/gpdd_results_v2.Rdata")

#increase threshold
predthreshold=0.2
gpdd_d2$predictable=select(gpdd_d2,R1a,R2a,R3a,R4a,R5a) %>% apply(1,FUN = function(x) ifelse(any(x>predthreshold), "yes","no"))
gpdd_d2$predictable_gr=select(gpdd_d2,R3m,R4m,R5m) %>% apply(1,FUN = function(x) ifelse(any(x>predthreshold), "yes","no"))
gpdd_d2$predictable_ag=ifelse(gpdd_d2$predictable=="yes" & gpdd_d2$predictable_gr=="yes", "ag",
                             ifelse(gpdd_d2$predictable=="yes" & gpdd_d2$predictable_gr=="no", "a",
                                    ifelse(gpdd_d2$predictable=="no" & gpdd_d2$predictable_gr=="yes", "g", "none")))
gpdd_d3=filter(gpdd_d2, predictable_ag!="none")
gpdd_results3=filter(gpdd_results, MainID %in% unique(gpdd_d3$MainID))



gpdd_d2$glebest= apply(gpdd_d2, 1, FUN=function(x) {unlist(x[paste0("gle",x$bestmodel)])})
gpdd_d2$lle_avgbest= apply(gpdd_d2, 1, FUN=function(x) {unlist(x[paste0("lle_avg",x$bestmodel)])})
gpdd_d2$lle_ppbest= apply(gpdd_d2, 1, FUN=function(x) {unlist(x[paste0("lle_pp",x$bestmodel)])})
gpdd_d2$glesign=ifelse(gpdd_d2$glebest>0, "positive", "negative")
gpdd_d2$lle_avgsign=ifelse(gpdd_d2$lle_avgbest>0, "positive", "negative")

focal_taxa=c("Aves","Osteichthyes", "Mammalia", "Bacillariophyceae", "Insecta")
gpdd_d2$TaxonomicClass2=ifelse(gpdd_d2$TaxonomicClass %in% focal_taxa, as.character(gpdd_d2$TaxonomicClass), "Other")

#temporary plots
plot(lle_avgbest~glebest, data=gpdd_d2, xlim=c(-1,1), ylim=c(-1,1))
abline(a = 0, b = 1)
abline(h=0)
abline(v=0)


#prop pos gle (0.33)
length(which(gpdd_d2$glebest>0))/length(which(!is.na(gpdd_d2$glebest)))
#gles estimated (154)
length(which(!is.na(gpdd_d2$glebest)))
#prop with any pos lle (0.48)
length(which(gpdd_d2$lle_ppbest>0))/length(which(!is.na(gpdd_d2$lle_ppbest)))
#prop with at least 0.5 pos lle (0.34)
length(which(gpdd_d2$lle_ppbest>0.5))/length(which(!is.na(gpdd_d2$lle_ppbest)))
#lle_avgs estimated (163)
length(which(!is.na(gpdd_d2$lle_avgbest)))
#prop pos avg lle (0.38)
length(which(gpdd_d2$lle_avgbest>0))/length(which(!is.na(gpdd_d2$lle_avgbest)))

plot(glebest~log(timescale_MinAge), data=gpdd_d2); abline(h=0)
plot(glebest~log(timescale_Lifespan), data=gpdd_d2); abline(h=0)
plot(glebest~log(Mass), data=gpdd_d2); abline(h=0)

plot(lle_ppbest~log(timescale1), data=gpdd_d2)
abline(h=0)

qplot(y = R1a, x= Reliability, data=gpdd_d2)
table(gpdd_d2$Reliability, gpdd_d2$predictability)
plot(gpdd_d2$Reliability, gpdd_d2$R1a, ylim=c(-1,1))
plot(gpdd_d$Reliability, gpdd_d$R1a, ylim=c(-1,1))

#counts of + and - by taxon, ag
signcounts=aggregate(glebest~glesign*TaxonomicClass2*predictable_ag, data=gpdd_d2, FUN=length) %>% 
  mutate(glesign2=ifelse(glesign=="negative", glebest*-1,glebest))
ggplot(signcounts, aes(x=TaxonomicClass2, y=glesign2)) + 
  facet_grid(predictable_ag~., scales = "free_y") + #geom_jitter(alpha=0.5, size=3, height = 0, width = 0.1) +
  geom_bar(aes(fill=glesign), stat = "identity") +
  geom_hline(yintercept = 0)  + theme_bw() + xlabvert + ylab("Count")

#histogram by ag, taxon
ggplot(gpdd_d3, aes(x=glebest, fill=TaxonomicClass2)) + 
  facet_grid(predictable_ag~.) + geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() + xlim(c(-2,2))

#by timescale, taxon, ag (don't have timescale data for all species)
ggplot(gpdd_d2, aes(y=glebest, x=log(timescale_MinAge), fill=TaxonomicClass2, color=glebest>0)) + 
  facet_grid(predictable_ag~.) + geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d2, aes(y=glebest, x=log(timestep_MinAge), fill=TaxonomicClass2, color=glebest>0)) + 
  facet_grid(predictable_ag~.) + geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d2, aes(y=log(timescale_MinAge), x=log(timestep_MinAge), fill=glesign)) + 
  facet_grid(predictable_ag~.) + geom_point(size=2.5, pch=21, color="black") +
  geom_vline(xintercept = 0) +
  classic 
ggplot(gpdd_d2, aes(y=log(timescale_Lifespan), x=log(timestep_Lifespan), fill=glesign)) + 
  facet_grid(predictable_ag~.) + geom_point(size=2.5, pch=21, color="black") +
  geom_vline(xintercept = 0) +
  classic 



ggplot(gpdd_d2, aes(x=lle_avgbest, fill=TaxonomicClass2)) + 
  facet_grid(predictable_ag~.) + geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw()

ggplot(gpdd_d2, aes(x=lle_ppbest)) + 
  facet_grid(predictable_ag~.) + geom_histogram(boundary = 0, binwidth = 0.05) + theme_bw()

test1=filter(gpdd_d2, glebest>2)
test=filter(gpdd_results2, gpdd_d2$gle1>5)
test1=filter(gpdd_d2, MainID==9921) #1188
test=filter(gpdd_results2, MainID==9921)
plot(test$modelresults1[[1]]$resultsdf$time, test$modelresults1[[1]]$resultsdf$obs, type="l")
lines(test$modelresults1[[1]]$resultsdf$time, test$modelresults1[[1]]$resultsdf$pred, col="red")
plot(test$modelresults3[[1]]$resultsdf$time, test$modelresults3[[1]]$resultsdf$obs, type="l")
lines(test$modelresults3[[1]]$resultsdf$time, test$modelresults3[[1]]$resultsdf$pred, col="red")
View(test1$data_rescale)
plot(test1$data_rescale[[1]]$SeriesStep, test1$data_rescale[[1]]$PopRescale_log, type="l")
lines(test$modelresults1[[1]]$resultsdf$time, log(test$modelresults1[[1]]$resultsdf$pred), col="red")

plot(test$modelresults2[[1]]$resultsdf$time, test$modelresults2[[1]]$resultsdf$obs, type="l")
lines(test$modelresults2[[1]]$resultsdf$time, test$modelresults2[[1]]$resultsdf$pred, col="red")

plot(test$modelresults2[[1]]$resultsdf$time, test$modelresults2[[1]]$resultsdf$Obs_abund, type="l")
lines(test$modelresults2[[1]]$resultsdf$time, test$modelresults2[[1]]$resultsdf$Pred_abund, col="red")

#Regression method

plot(gpdd_d2$gle4, gpdd_d2$gle5)
abline(a = 0, b=1)
abline(h=0)
abline(v=0)
