# GPDD stability analysis - plots

library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")
library("purrr")

load(file = "./data/gpdd_results_fixedEtau.Rdata")
source("~/GRAD SCHOOL/R reference/ggplot themes rogers.R")

#general stats

#number of distinct taxa (140)
n_distinct(gpdd_d$TaxonID)
#gles estimated (171)
length(which(!is.na(gpdd_d$gle)))
#prop pos gles 
length(which(gpdd_d$gle>0.01))/length(which(!is.na(gpdd_d$gle)))
#prop pos gles 1d 
length(which(gpdd_d$gle1d>0.01))/length(which(!is.na(gpdd_d$gle1d)))
#prop pos LEshift
length(which(gpdd_d$minci>0.01))/length(which(!is.na(gpdd_d$minci)))
#prop pos LEshift 1d
length(which(gpdd_d$minci1d>0.01))/length(which(!is.na(gpdd_d$minci1d)))

#prop with any pos lle 
length(which(gpdd_d$lle_pp>0))/length(which(!is.na(gpdd_d$lle_pp)))
#prop with at least 0.5 pos lle 
length(which(gpdd_d$lle_pp>0.5))/length(which(!is.na(gpdd_d$lle_pp)))

# #gles are correlated, although some outliers
# pairs(select(gpdd_d, gle1:gle5, gle))
# #timespan variable are correlated with each other, or their inverse
# pairs(select(gpdd_d, MinAge_mo, Lifespan_mo, timescale_MinAge, timescale_Lifespan, timestep_MinAge, timestep_Lifespan))

#gle is never higher than lle_avg
plot(lle_avgbest~gle, data=gpdd_d, xlim=c(-1,1), ylim=c(-1,1))
abline(a = 0, b = 1); abline(h=0); abline(v=0)
#gle vs minci (mostly correlated)
plot(minci~gle, data=gpdd_d, xlim=c(-1,1), ylim=c(-1,1))
abline(a = 0, b = 1); abline(h=0); abline(v=0)

#GLE vs. LLE prop positive
ggplot(gpdd_d, aes(x=lle_pp, fill=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.05) + theme_bw()
ggplot(gpdd_d, aes(x=lle_pp, fill=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.05) + theme_bw()

#reliability seems to be unrelated to predictability
table(gpdd_d$Reliability, gpdd_d$predictable_ag)
plot(gpdd_d$Reliability, gpdd_d$bestR2, ylim=c(-1,1))

#counts of gle signs by taxon (predictability)
signcounts=aggregate(minci~LEshiftsign*TaxonomicClass2*predictable_ag, data=gpdd_d, FUN=length) %>% 
  mutate(LEshiftsign2=ifelse(LEshiftsign=="not chaotic", minci*-1,minci))
ggplot(signcounts, aes(x=TaxonomicClass2, y=LEshiftsign2)) + 
  #facet_grid(predictable_ag~., scales = "free_y") + #geom_jitter(alpha=0.5, size=3, height = 0, width = 0.1) +
  geom_bar(aes(fill=LEshiftsign), stat = "identity") +
  geom_hline(yintercept = 0)  + theme_bw() + xlabvert + ylab("Count") + ylim(c(-55,20))
#compare to 1d
signcounts1d=aggregate(minci1d~LEshiftsign1d*TaxonomicClass2*predictable_ag, data=gpdd_d, FUN=length) %>% 
  mutate(LEshiftsign1d2=ifelse(LEshiftsign1d=="not chaotic", minci1d*-1,minci1d))
ggplot(signcounts1d, aes(x=TaxonomicClass2, y=LEshiftsign1d2)) + 
  #facet_grid(predictable_ag~., scales = "free_y") + #geom_jitter(alpha=0.5, size=3, height = 0, width = 0.1) +
  geom_bar(aes(fill=LEshiftsign1d), stat = "identity") +
  geom_hline(yintercept = 0)  + theme_bw() + xlabvert + ylab("Count") + ylim(c(-55,20))

#gle vs gle1d
ggplot(filter(gpdd_d, E>1), aes(y=minci1d, x=minci, color=E)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2,stroke=1.5) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  classic 

#histogram gle by taxon (predictability)
ggplot(gpdd_d, aes(x=minci, fill=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() 
ggplot(gpdd_d, aes(x=minci1d, fill=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() 
ggplot(gpdd_d, aes(x=les_mo, fill=TaxonomicClass2)) + 
  facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() 
ggplot(gpdd_d, aes(x=les_gen, fill=TaxonomicClass2)) + 
  facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() + xlim(c(-2,2))
#by sampling interval
ggplot(gpdd_d, aes(x=minci, fill=TaxonomicClass2)) + 
  facet_grid(SamplingInterval~.) + 
  geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() 

#distribution of Es, taus
ggplot(gpdd_d, aes(x=factor(E), fill=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  geom_bar(position = "stack") + xlab("Best E") +
  classic
ggplot(gpdd_d, aes(y=E, x=datasetlength, color=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point() + 
  classic
ggplot(gpdd_d, aes(x=factor(tau), fill=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  geom_bar(position = "stack") + xlab("Best tau") +
  classic

#gle by E, theta
ggplot(gpdd_d, aes(y=minci, x=factor(E), fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) + xlab("Best E") +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=minci, x=theta, fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 

#gle by latitude
ggplot(gpdd_d, aes(y=minci, x=abs(LatDD), fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  #facet_wrap(~TaxonomicClass2, nrow=2) +
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=les_mo, x=abs(LatDD), fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  #facet_wrap(~TaxonomicClass2, nrow=2) +
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=les_gen, x=abs(LatDD), fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  #facet_wrap(~TaxonomicClass2, nrow=2) +
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic + ylim(c(-3,3))

#gle by intrinsic timescale
ggplot(gpdd_d, aes(y=les_mo, x=log10(MinAge_mo), fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=les_gen, x=log10(MinAge_mo), fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic + ylim(c(-3,3))
ggplot(gpdd_d, aes(y=les_mo, x=log10(Lifespan_mo), fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 

#gle by timescale ratio
#timescale_MinAge = gens/ts length
#timestep_MinAge = gens/timestep
ggplot(gpdd_d, aes(y=minci, x=log10(timestep_MinAge), fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) + xlab("log10 gens/timestep") +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=les_mo, x=log10(timestep_MinAge), fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) + xlab("log10 gens/timestep") +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=les_gen, x=log10(timestep_MinAge), fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) + xlab("log10 gens/timestep") +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic + ylim(c(-3,3))
ggplot(gpdd_d, aes(y=gle_mo, x=log10(timescale_MinAge), fill=TaxonomicClass2, color=LEshiftsign)) + 
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
ggplot(gpdd_d, aes(y=bestR2, x=log10(timestep_MinAge), fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) + xlab("log10 gens/timestep") +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=monotonicR2, x=log10(timestep_MinAge), fill=TaxonomicClass2, color=LEshiftsign)) + 
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
ggplot(gpdd_d, aes(y=minci, x=monotonicR2, fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=les_mo, x=monotonicR2, fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 

#gle vs mass
ggplot(gpdd_d, aes(y=les_mo, x=log10(Mass_g), fill=TaxonomicClass2, color=LEshiftsign)) + 
  #facet_grid(predictable_ag~.) +
  #facet_grid(TaxonomicClass2~., scales="free_y") +
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("black", "red", NA)) +
  classic 
ggplot(gpdd_d, aes(y=gle_gen, x=log10(Mass_g), fill=TaxonomicClass2, color=LEshiftsign)) + 
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
high=filter(gpdd_d, gle>1)
par(mfrow=c(2,3))
for(i in 1:nrow(high)) {
  dtemp=high$data_rescale[i][[1]]
  plot(PopRescale_log~SeriesStep, data=dtemp, type="l", main=paste(i, high$CommonName[i]))
}

#examine duplicates
duplicates=aggregate(gle~TaxonID, data=gpdd_d, FUN=length) %>% filter(gle>1)
dup=filter(gpdd_d, TaxonID %in% duplicates$TaxonID) %>% arrange(TaxonID)
dupsimp=select(dup, MainID, TaxonID, CommonName, ExactName, Country, SamplingUnits, gle, glesign)
# Of 21 species with multiple ts, 7 species have LEs of different signs
# includes voles in diff areas, dramatically different LEs
# write.csv(dupsimp, "./data/duplicatespecies.csv", row.names = F)

missing=filter(gpdd_d, is.na(gle))
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

ser=hightau[1,]$data_rescale[[1]]$PopRescale_log