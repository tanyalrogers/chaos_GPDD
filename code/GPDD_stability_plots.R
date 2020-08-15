# GPDD stability analysis - plots

library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")
library("purrr")

load(file = "./data/gpdd_results.Rdata")
source("~/GRAD SCHOOL/R reference/ggplot themes rogers.R")

#general stats

#number of distinct taxa (140)
n_distinct(gpdd_d$TaxonID)
#gles estimated (169)
length(which(!is.na(gpdd_d$glebest)))
#prop pos gles (0.32)
length(which(gpdd_d$glebest>0))/length(which(!is.na(gpdd_d$glebest)))
#prop pos gles 1d (0.19)
length(which(gpdd_d$glebest1d>0))/length(which(!is.na(gpdd_d$glebest1d)))

#lle_avgs estimated (177)
length(which(!is.na(gpdd_d$lle_avgbest)))
#prop pos avg lle (0.37)
length(which(gpdd_d$lle_avgbest>0))/length(which(!is.na(gpdd_d$lle_avgbest)))

#prop with any pos lle (0.71)
length(which(gpdd_d$lle_ppbest>0))/length(which(!is.na(gpdd_d$lle_ppbest)))
#prop with at least 0.5 pos lle (0.39)
length(which(gpdd_d$lle_ppbest>0.5))/length(which(!is.na(gpdd_d$lle_ppbest)))

#gles are correlated, although some outliers
pairs(select(gpdd_d, gle1:gle5, glebest))
pairs(select(gpdd_d, MinAge_mo, Lifespan_mo, timescale_MinAge, timescale_Lifespan, timestep_MinAge, timestep_Lifespan))

#gle is never higher than lle_avg
plot(lle_avgbest~glebest, data=gpdd_d, xlim=c(-1,1), ylim=c(-1,1))
abline(a = 0, b = 1); abline(h=0); abline(v=0)

#LLE prop positive
ggplot(gpdd_d, aes(x=lle_ppbest, fill=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.05) + theme_bw()

#reliability seems to be unrelated to predictability
table(gpdd_d$Reliability, gpdd_d$predictable_ag)
plot(gpdd_d$Reliability, gpdd_d$R1a, ylim=c(-1,1))
plot(gpdd_d$Reliability, gpdd_d$R1a, ylim=c(-1,1))

#counts of gle signs by taxon (predictability)
signcounts=aggregate(glebest~glesign*TaxonomicClass2*predictable_ag, data=gpdd_d, FUN=length) %>% 
  mutate(glesign2=ifelse(glesign=="negative", glebest*-1,glebest))
ggplot(signcounts, aes(x=TaxonomicClass2, y=glesign2)) + 
  #facet_grid(predictable_ag~., scales = "free_y") + #geom_jitter(alpha=0.5, size=3, height = 0, width = 0.1) +
  geom_bar(aes(fill=glesign), stat = "identity") +
  geom_hline(yintercept = 0)  + theme_bw() + xlabvert + ylab("Count") + ylim(c(-55,20))
#compare to 1d
signcounts1d=aggregate(glebest1d~glesign1d*TaxonomicClass2*predictable_ag, data=gpdd_d, FUN=length) %>% 
  mutate(glesign2=ifelse(glesign1d=="negative", glebest1d*-1,glebest1d))
ggplot(signcounts1d, aes(x=TaxonomicClass2, y=glesign2)) + 
  #facet_grid(predictable_ag~., scales = "free_y") + #geom_jitter(alpha=0.5, size=3, height = 0, width = 0.1) +
  geom_bar(aes(fill=glesign1d), stat = "identity") +
  geom_hline(yintercept = 0)  + theme_bw() + xlabvert + ylab("Count") + ylim(c(-55,20))

#gle vs gle1d
ggplot(filter(gpdd_d, Ebest>1), aes(y=glebest1d, x=glebest, color=Ebest)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2,stroke=1.5) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  classic 

#histogram gle by taxon (predictability)
ggplot(gpdd_d, aes(x=glebest, fill=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() 
ggplot(gpdd_d, aes(x=glebest1d, fill=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() 
ggplot(gpdd_d, aes(x=gle_mo, fill=TaxonomicClass2)) + 
  facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() 
ggplot(gpdd_d, aes(x=gle1d_mo, fill=TaxonomicClass2)) + 
  facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() + xlim(c(-2,2))
ggplot(gpdd_d, aes(x=gle_gen, fill=TaxonomicClass2)) + 
  facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() + xlim(c(-2,2))
#by sampling interval
ggplot(gpdd_d, aes(x=glebest, fill=TaxonomicClass2)) + 
  facet_grid(SamplingInterval~.) + 
  geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() 

#distribution of Es
ggplot(gpdd_d, aes(x=factor(Ebest), fill=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  geom_bar(position = "stack") + xlab("Best E") +
  classic
ggplot(gpdd_d, aes(y=Ebest, x=datasetlength, color=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point() + 
  classic

#gle by E, theta
ggplot(gpdd_d, aes(y=glebest, x=factor(Ebest), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) + xlab("Best E") +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=glebest, x=thetabest, fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 

#gle by latitude
ggplot(gpdd_d, aes(y=glebest, x=abs(LatDD), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  #facet_wrap(~TaxonomicClass2, nrow=2) +
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=gle_mo, x=abs(LatDD), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  #facet_wrap(~TaxonomicClass2, nrow=2) +
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=gle_gen, x=abs(LatDD), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  #facet_wrap(~TaxonomicClass2, nrow=2) +
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic + ylim(c(-3,3))

#gle by intrinsic timescale
ggplot(gpdd_d, aes(y=gle_mo, x=log10(MinAge_mo), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=gle_gen, x=log10(MinAge_mo), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic + ylim(c(-3,3))
ggplot(gpdd_d, aes(y=gle_mo, x=log10(Lifespan_mo), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 

#gle by timescale ratio
#timescale_MinAge = gens/ts length
#timestep_MinAge = gens/timestep
ggplot(gpdd_d, aes(y=glebest, x=log10(timestep_MinAge), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) + xlab("log10 gens/timestep") +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=gle_mo, x=log10(timestep_MinAge), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) + xlab("log10 gens/timestep") +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=gle_gen, x=log10(timestep_MinAge), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) + xlab("log10 gens/timestep") +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic + ylim(c(-3,3))
ggplot(gpdd_d, aes(y=gle_mo, x=log10(timescale_MinAge), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) +
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=gle_gen, x=log10(timescale_MinAge), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) +
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic + 
  ylim(c(-3,3))

#R2 by gen/timestep
ggplot(gpdd_d, aes(y=R2abundbest, x=log10(timestep_MinAge), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) + xlab("log10 gens/timestep") +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=monotonicR2, x=log10(timestep_MinAge), fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) + xlab("log10 gens/timestep") +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 

#gens/timespan vs gens/timestep
ggplot(gpdd_d, aes(y=log10(timescale_MinAge), x=log10(timestep_MinAge), fill=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2.5, pch=21, color="black") +
  geom_vline(xintercept = 0) +
  classic 
ggplot(gpdd_d, aes(y=log10(timescale_Lifespan), x=log10(timestep_Lifespan), fill=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2.5, pch=21, color="black") +
  geom_vline(xintercept = 0) +
  classic 

#gle vs monotonic trend
ggplot(gpdd_d, aes(y=glebest, x=monotonicR2, fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=gle_mo, x=monotonicR2, fill=TaxonomicClass2, color=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 

#gle vs mass
ggplot(gpdd_d, aes(y=gle_mo, x=log10(Mass_g), fill=TaxonomicClass2, color=glebest>0)) + 
  #facet_grid(predictable_ag~.) +
  #facet_grid(TaxonomicClass2~., scales="free_y") +
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=gle_gen, x=log10(Mass_g), fill=TaxonomicClass2, color=glebest>0)) + 
  #facet_grid(predictable_ag~.) +
  facet_grid(TaxonomicClass2~., scales="free_y") +
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic #+ ylim(c(-3,3))

#diagnostics

test1=filter(gpdd_d, MainID==9921)
test=filter(gpdd_results, MainID==9921)
#plot obs and pred
plot(test$modelresultsbest[[1]]$resultsdf$time, test$modelresultsbest[[1]]$resultsdf$obs, type="l")
lines(test$modelresultsbest[[1]]$resultsdf$time, test$modelresultsbest[[1]]$resultsdf$pred, col="red")
#plot obs and pred abund
plot(test$modelresultsbest[[1]]$resultsdf$time, test$modelresultsbest[[1]]$resultsdf$Obs_abund, type="l")
lines(test$modelresultsbest[[1]]$resultsdf$time, test$modelresultsbest[[1]]$resultsdf$Pred_abund, col="red")

#high values
high=filter(gpdd_d, glebest>1)
par(mfrow=c(2,3))
for(i in 1:nrow(high)) {
  dtemp=high$data_rescale[i][[1]]
  plot(PopRescale_log~SeriesStep, data=dtemp, type="l", main=paste(i, high$CommonName[i]))
}

#examine duplicates
duplicates=aggregate(glebest~TaxonID, data=gpdd_d, FUN=length) %>% filter(glebest>1)
dup=filter(gpdd_d, TaxonID %in% duplicates$TaxonID) %>% arrange(TaxonID)
dupsimp=select(dup, MainID, TaxonID, CommonName, ExactName, Country, SamplingUnits, glebest, glesign)
# Of 21 species with multiple ts, 7 species have LEs of different signs
# includes voles in diff areas, dramatically different LEs
# write.csv(dupsimp, "./data/duplicatespecies.csv", row.names = F)

missing=filter(gpdd_d, is.na(glebest))
fish=filter(gpdd_d, TaxonomicClass2=="Osteichthyes")
diatoms=filter(gpdd_d, TaxonomicClass2=="Bacillariophyceae")

for(i in 1:nrow(diatoms)) {
  dtemp=diatoms$data_rescale[i][[1]]
  plot(PopRescale_log~SeriesStep, data=dtemp, type="l", main=paste(i, diatoms$CommonName[i]))
}

#for plotting 
for(i in 1:8) {
  dtemp=missing$data_rescale[i][[1]]
  plot(PopRescale_log~SeriesStep, data=dtemp, type="l", main=paste(i, missing$CommonName[i]))
}

cor(dup$data_rescale[53][[1]]$PopRescale_log,dup$data_rescale[54][[1]]$PopRescale_log)
