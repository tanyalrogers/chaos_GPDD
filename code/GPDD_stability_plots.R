# GPDD stability analysis - plots

library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")
library("purrr")

load(file = "./data/gpdd_results_update.Rdata")
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
#prop pos minci
length(which(gpdd_d$minci>0.01))/length(which(!is.na(gpdd_d$minci)))
#prop pos minci 1d
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
# plot(lle_avgbest~gle, data=gpdd_d, xlim=c(-1,1), ylim=c(-1,1))
# abline(a = 0, b = 1); abline(h=0); abline(v=0)
#gle vs minci (mostly correlated)
plot(minci~gle, data=gpdd_d, xlim=c(-1,1), ylim=c(-1,1))
abline(a = 0, b = 1); abline(h=0); abline(v=0)

#GLE vs. LLE prop positive
ggplot(gpdd_d, aes(x=lle_pp, fill=glesign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.05) + classic
ggplot(gpdd_d, aes(x=lle_pp, fill=mincisign)) + 
  #facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.05) + classic

#reliability seems to be unrelated to predictability
table(gpdd_d$Reliability, gpdd_d$predictable_ag)
plot(gpdd_d$Reliability, gpdd_d$bestR2, ylim=c(-1,1))

#classifications by taxon (predictability)
signcountsnd=aggregate(minci~mincisign*TaxonomicClass2*predictable_ag, data=gpdd_d, FUN=length) %>% 
  mutate(mincisign2=ifelse(mincisign=="not chaotic", minci*-1,minci), Econ="Free E")
signcounts1d=aggregate(minci1d~mincisign1d*TaxonomicClass2*predictable_ag, data=gpdd_d, FUN=length) %>% 
  mutate(mincisign1d2=ifelse(mincisign1d=="not chaotic", minci1d*-1,minci1d), Econ="E = 1")
colnames(signcounts1d)<-colnames(signcounts)
signcounts=rbind(signcountsnd,signcounts1d); signcounts$Econ=factor(signcounts$Econ, levels=unique(signcounts$Econ))
ggplot(signcounts, aes(x=TaxonomicClass2, y=mincisign2)) + 
  #facet_grid(predictable_ag~., scales = "free_y") + #geom_jitter(alpha=0.5, size=3, height = 0, width = 0.1) +
  facet_grid(.~Econ) + 
  geom_bar(aes(fill=mincisign), stat = "identity") +
  geom_hline(yintercept = 0)  + classic + xlabvert + ylab("Count") +
  labs(fill="Classification", x="Taxonomic Class") + removefacetbackground
# #compare to 1d
# ggplot(signcounts1d, aes(x=TaxonomicClass2, y=mincisign1d2)) + 
#   #facet_grid(predictable_ag~., scales = "free_y") + #geom_jitter(alpha=0.5, size=3, height = 0, width = 0.1) +
#   geom_bar(aes(fill=mincisign1d), stat = "identity") +
#   geom_hline(yintercept = 0)  + classic + xlabvert + ylab("Count") + ylim(c(-60,20)) +
#   labs(fill="Classification", x="Taxonomic Class", title="E = 1")

#classifications by E
signcountsE=aggregate(minci~mincisign*E*TaxonomicClass2, data=gpdd_d, FUN=length) %>% 
  mutate(mincisign2=ifelse(mincisign=="not chaotic", minci*-1,minci)) %>% 
  group_by(E) %>% mutate(signprop=minci/sum(minci))
ggplot(signcountsE, aes(x=factor(E), y=mincisign2)) + 
  facet_grid(TaxonomicClass2~., scales = "free_y") +
  geom_bar(aes(fill=mincisign), stat = "identity") +
  geom_hline(yintercept = 0)  + classic + ylab("Count") +
  labs(fill="Classification", x="E")
ggplot(signcountsE, aes(x=factor(E), y=minci)) + 
  facet_grid(TaxonomicClass2~., scales = "free_y") +
  geom_bar(aes(fill=mincisign), stat = "identity") +
  geom_hline(yintercept = 0)  + classic + ylab("Count") +
  labs(fill="Classification", x="E")
ggplot(signcountsE, aes(x=factor(E), y=signprop)) + 
  geom_bar(aes(fill=mincisign), stat = "identity") +
  classic + ylab("Proportion") + scale_y_continuous(expand = expand_scale(mult = c(0, 0))) +
  labs(fill="Classification", x="E") 

#gle vs gle1d
ggplot(filter(gpdd_d, E>1), aes(y=minci1d, x=minci, color=E)) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=2,stroke=1.5) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  classic + labs(x="LE lower bound - Free E", y="LE lower bound - E = 1", color="Free E")

#histogram gle by taxon (predictability)
ggplot(gpdd_d, aes(x=minci, fill=TaxonomicClass2)) + 
  facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0.01, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() + 
  labs(x="LE lower bound", y="Count", fill="Taxonomic\nClass") 
ggplot(gpdd_d, aes(x=minci1d, fill=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() 
ggplot(gpdd_d, aes(x=minci_mo, fill=TaxonomicClass2)) + 
  facet_grid(predictable_ag~.) + 
  geom_histogram(boundary = 0, binwidth = 0.1) +
  geom_vline(xintercept = 0) + theme_bw() 
ggplot(gpdd_d, aes(x=minci_gen, fill=TaxonomicClass2)) + 
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
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Taxonomic\nClass") 
ggplot(gpdd_d, aes(x=factor(tau), fill=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  geom_bar(position = "stack") + xlab("Best tau") +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Taxonomic\nClass")
ggplot(gpdd_d, aes(x=factor(tau), fill=factor(E))) + 
  #facet_grid(predictable_ag~.) + 
  geom_bar(position = "stack") + xlab("Best tau") +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Best E")

#E vs time series length
ggplot(gpdd_d, aes(y=E, x=datasetlength, color=TaxonomicClass2)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("E") + xlab("Time Series Length (timesteps)") + 
  geom_point(size=2, alpha=0.4) + 
  classic + labs(color="Taxonomic\nClass") + legalpha
ggplot(gpdd_d, aes(y=E, x=timescale_MinAge*MinAge_mo/12, color=TaxonomicClass2)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("E") + xlab("Time Series Length (years)") + 
  geom_point(size=2, alpha=0.4) + 
  classic + labs(color="Taxonomic\nClass") + legalpha
ggplot(gpdd_d, aes(y=E, x=log10(timescale_MinAge), color=TaxonomicClass2)) + 
  facet_grid(TaxonomicClass2~.) + 
  ylab("E") + xlab("Time Series Length (log10 generations)") + 
  geom_point(size=2, alpha=0.4) + 
  classic + labs(color="Taxonomic\nClass") + legalpha

#gle by E, theta
ggplot(gpdd_d, aes(y=minci, x=E)) + 
  xlab("E") + ylab("LE lower bound") +
  geom_point(aes(color=TaxonomicClass2), size=2, alpha=0.5) + 
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  geom_quantile(quantiles = c(0.20, 0.5, 0.80), formula=y~poly(x,2), size=2, alpha=0.5) +
  classic + labs(color="Taxonomic\nClass") + legalpha
ggplot(gpdd_d, aes(y=minci, x=factor(E), color=TaxonomicClass2)) + 
  #facet_grid(TaxonomicClass2~., scales="free_y") + 
  xlab("E") + ylab("LE lower bound") +
  geom_point(size=2, alpha=0.5) + 
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass") + legalpha
ggplot(gpdd_d, aes(y=minci, x=theta, color=TaxonomicClass2)) + 
  xlab("Theta") + ylab("LE lower bound") +
  geom_point(size=2, alpha=0.5) + 
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass") + legalpha

#age at maturity vs lifespan and mass
ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=log10(Lifespan_mo), fill=datasetlength)) + 
  ylab("log10 Age at Maturity (months)") + xlab("log10 Lifespan (months)") + 
  geom_point(size=2.5, pch=21, color="black") +
  classic + labs(fill="Time Series\nLength")
ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=log10(Mass_g), fill=datasetlength)) + 
  ylab("log10 Age at Maturity (months)") + xlab("log10 Mass (g)") +
  geom_point(size=2.5, pch=21, color="black") +
  classic + labs(fill="Time Series\nLength")
ggplot(gpdd_d, aes(y=log10(Lifespan_mo), x=log10(Mass_g), fill=datasetlength)) + 
  ylab("log10 Lifespan (months)") + xlab("log10 Mass (g)") +
  geom_point(size=2.5, pch=21, color="black") +
  classic + labs(fill="Time Series\nLength")
#gens/timespan vs gens/timestep
ggplot(gpdd_d, aes(y=log10(timestep_MinAge), x=log10(MinAge_mo), fill=datasetlength)) + 
  ylab("log10 Generations/Timestep") + xlab("log10 Age at Maturity (months)") +
  geom_point(size=2.5, pch=21, color="black") +
  classic + labs(fill="Time Series\nLength")
ggplot(gpdd_d, aes(y=log10(timescale_MinAge), x=log10(timestep_MinAge), fill=datasetlength)) + 
  ylab("log10 Generations sampled") + xlab("log10 Generations/Timestep") +
  geom_point(size=2.5, pch=21, color="black") +
  classic + labs(fill="Time Series\nLength")
ggplot(gpdd_d, aes(y=log10(timescale_Lifespan), x=log10(timestep_Lifespan), fill=datasetlength)) + 
  ylab("log10 Lifespans sampled") + xlab("log10 Lifespans/Timestep") +
  geom_point(size=2.5, pch=21, color="black") +
  classic + labs(fill="Time Series\nLength")

#gen time vs. ts length
ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=datasetlength, color=TaxonomicClass2)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("log10 Generation time (months)") + xlab("Time Series Length (timesteps)") + 
  geom_point(size=2, alpha=0.4) + 
  classic + labs(color="Taxonomic\nClass") + legalpha
ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=timescale_MinAge*MinAge_mo/12, color=TaxonomicClass2)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("log10 Generation time (months)") + xlab("Time Series Length (years)") + 
  geom_point(size=2, alpha=0.4) + 
  classic + labs(color="Taxonomic\nClass") + legalpha
ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=log10(timescale_MinAge), color=TaxonomicClass2)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("log10 Generation time (months)") + xlab("Time Series Length (log10 generations)") + 
  geom_point(size=2, alpha=0.4) + 
  classic + labs(color="Taxonomic\nClass") + legalpha

#gle vs mass
ggplot(gpdd_d, aes(y=minci_mo, x=log10(Mass_g), color=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) +
  #facet_grid(TaxonomicClass2~., scales="free_y") +
  ylab("LE lower bound (per month)") + xlab("log10 Mass (g)") +
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")
#comparison to Anderson and Gilooly 2020 (positive values only)
ggplot(filter(gpdd_d, mincisign=="chaotic"), aes(y=log10(minci_mo), x=log10(Mass_g))) + 
  #facet_grid(TaxonomicClass2~., scales="free_y") +
  ylab("log10 LE lower bound (per month)") + xlab("log10 Mass (g)") +
  geom_smooth(method="lm", se = F, color="black") +
  geom_point(aes(color=TaxonomicClass2), size=2, alpha=0.5) + 
  classic + labs(color="Taxonomic\nClass") + legalpha
summary(lm(log10(minci_mo)~log10(Mass_g), data=filter(gpdd_d, mincisign=="chaotic")))
ggplot(filter(gpdd_d, mincisign=="chaotic"), aes(y=log10(minci_gen), x=log10(Mass_g))) + 
  #facet_grid(TaxonomicClass2~., scales="free_y") +
  ylab("log10 LE lower bound (per generation)") + xlab("log10 Mass (g)") +
  geom_smooth(method="lm", se = F, color="black") +
  geom_hline(yintercept = log10(0.06), lty=2) +
  geom_point(aes(color=TaxonomicClass2), size=2, alpha=0.5) + 
  classic + labs(color="Taxonomic\nClass") + legalpha
summary(lm(log10(minci_gen)~log10(Mass_g), data=filter(gpdd_d, mincisign=="chaotic")))

ggplot(gpdd_d, aes(y=minci_gen, x=log10(Mass_g), color=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) +
  #facet_grid(TaxonomicClass2~., scales="free_y") +
  ylab("LE lower bound (per generation)") + xlab("log10 Mass (g)") +
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass") + ylim(c(-.3,.3)) #ylim(c(-3,3))

#gle by gen time
ggplot(gpdd_d, aes(y=minci_mo, x=log10(MinAge_mo), color=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  ylab("LE lower bound (per month)") + xlab("log10 Generation Time (months)") + 
  geom_point(size=2, alpha=0.4) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass") + legalpha
ggplot(gpdd_d, aes(y=minci_gen, x=log10(MinAge_mo), color=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  ylab("LE lower bound (per generation)") + xlab("log10 Generation Time (months)") + 
  geom_point(size=2, alpha=0.4) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass") + ylim(c(-.3,.3)) + xlim(c(-2,0))
#gle by gen time by E
ggplot(gpdd_d, aes(y=minci_mo, x=log10(MinAge_mo), color=E)) + 
  facet_grid(TaxonomicClass2~., scales="free_y") + 
  ylab("LE lower bound (per month)") + xlab("log10 Generation Time (months)") + 
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic #+ labs(color="Taxonomic\nClass")

#gle by ts length
ggplot(gpdd_d, aes(y=minci_mo, x=datasetlength, color=TaxonomicClass2)) + 
  #facet_grid(TaxonomicClass2~., scales="free_y") + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (timesteps)") +
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")
ggplot(gpdd_d, aes(y=minci_mo, x=log10(timescale_MinAge), color=TaxonomicClass2)) + 
  #facet_grid(TaxonomicClass2~., scales="free_y") + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")
#gle by ts length by E
ggplot(gpdd_d, aes(y=minci_mo, x=log10(timescale_MinAge), color=E)) + 
  #facet_grid(TaxonomicClass2~., scales="free_y") + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic #+ labs(color="Taxonomic\nClass")
#gle by ts length by gen time
ggplot(gpdd_d, aes(y=minci_mo, x=log10(timescale_MinAge), color=log10(MinAge_mo))) + 
  #facet_grid(TaxonomicClass2~., scales="free_y") + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="log10\nGeneration\ntime (months)")

#R2 by time series length
ggplot(gpdd_d, aes(y=bestR2, x=log10(timescale_MinAge), color=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  ylab("R-squared") + xlab("Time Series Length (log10 Generations)") +
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass") + legalpha

#gle by R2
ggplot(gpdd_d, aes(y=minci, x=bestR2, color=TaxonomicClass2)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("LE lower bound (per timestep)") + xlab("R-squared") +
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")
ggplot(gpdd_d, aes(y=minci, x=bestR2m, color=TaxonomicClass2)) + 
  facet_grid(TaxonomicClass2~.) + 
  ylab("LE lower bound (per timestep)") + xlab("R-squared (growth rate)") +
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")

#gle vs monotonic trend
ggplot(gpdd_d, aes(y=minci, x=monotonicR2, color=TaxonomicClass2)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("LE lower bound (per timestep)") + xlab("Monotonic trend R-squared") +
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Classification", fill="Taxonomic\nClass")
ggplot(gpdd_d, aes(y=minci_mo, x=monotonicR2, color=TaxonomicClass2)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("LE lower bound (per month)") + xlab("Monotonic trend R-squared") +
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Classification", fill="Taxonomic\nClass")
#ts length vs monotonic trend
ggplot(gpdd_d, aes(y=monotonicR2, x=log10(timestep_MinAge), color=TaxonomicClass2)) + 
  #facet_grid(predictable_ag~.) + 
  ylab("Monotonic trend R-squared") + xlab("Time Series Length (log10 Generations)") +
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")
#R2 vs monotonic trend
ggplot(gpdd_d, aes(y=bestR2, x=monotonicR2, color=E)) + 
  facet_grid(TaxonomicClass2~.) + 
  ylab("R-squared") + xlab("Monotonic trend R-squared") +
  geom_point(size=2, alpha=0.5) +
  geom_abline(slope = 1, intercept = 0) + scale_color_viridis_c() +
  classic #+ labs(color="Taxonomic\nClass")
#E vs monotonic trend
ggplot(gpdd_d, aes(y=E, x=monotonicR2, color=mincisign)) + 
  facet_grid(TaxonomicClass2~.) + 
  ylab("E") + xlab("Monotonic trend R-squared") +
  geom_point(size=2, alpha=0.5) + #scale_color_viridis_c() +
  classic + labs(color="Classification")

#gle vs latitude
ggplot(gpdd_d, aes(y=minci, x=abs(LatDD), fill=TaxonomicClass2, color=mincisign)) + 
  facet_wrap(~TaxonomicClass2, nrow=2) +
  ylab("LE lower bound") + xlab("Latitude (abs. value)") + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Classification", fill="Taxonomic\nClass")
ggplot(gpdd_d, aes(y=minci_mo, x=abs(LatDD), fill=TaxonomicClass2, color=mincisign)) + 
  #facet_wrap(~TaxonomicClass2, nrow=2) +
  ylab("LE lower bound (per month)") + xlab("Latitude (abs. value)") + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Classification", fill="Taxonomic\nClass")
ggplot(gpdd_d, aes(y=minci_gen, x=abs(LatDD), fill=TaxonomicClass2, color=mincisign)) + 
  #facet_wrap(~TaxonomicClass2, nrow=2) +
  ylab("LE lower bound (per generation)") + xlab("Latitude (abs. value)") +
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Classification", fill="Taxonomic\nClass") +
  ylim(c(-3,3))

#gle vs growth rate
gpdd_d$meangr=map_dbl(gpdd_d$data_rescale, ~mean(.x$PopRescale_gr, na.rm=T))
gpdd_d$maxgr=map_dbl(gpdd_d$data_rescale, ~max(.x$PopRescale_gr, na.rm=T))
gpdd_d$maxgr_mo=gpdd_d$maxgr/timescale_mo(gpdd_d$SamplingInterval,1)

ggplot(gpdd_d, aes(y=minci_mo, x=meangr, color=TaxonomicClass2)) + 
  facet_grid(.~E) +
  ylab("LE lower bound (per month)") + xlab("Mean growth rate") + 
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")
ggplot(gpdd_d, aes(y=minci_mo, x=maxgr, color=TaxonomicClass2)) + 
  #facet_grid(.~E) +
  ylab("LE lower bound (per month)") + xlab("Max growth rate") + 
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")

ggplot(gpdd_d, aes(y=log10(maxgr_mo), x=log10(Mass_g), color=mincisign)) + 
  #facet_grid(.~E) +
  ylab("log10 Max growth rate (monthly)") + xlab("log10 Mass(g)") + 
  geom_point(size=2, alpha=0.5) +
  classic + labs(color="Classification")


#resolving tslength vs gen time
gpdd_d$log10_MinAge_mo_bin=cut(log10(gpdd_d$MinAge_mo), seq(from=-2, to=2.5, by=0.5))
ggplot(gpdd_d, aes(y=minci_mo, x=log10(timescale_MinAge), color=TaxonomicClass2)) + 
  facet_wrap(log10_MinAge_mo_bin~., nrow = 3, scales="free") + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")
ggplot(filter(gpdd_d, MinAge_mo %in% c(8,10,12,24)), aes(y=minci_mo, x=timescale_MinAge*MinAge_mo/12, color=TaxonomicClass2)) + 
  facet_grid(MinAge_mo~., scales="free") + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (years)") + 
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")
ggplot(filter(gpdd_d, MinAge_mo %in% c(8,10,12,24)), aes(y=minci_mo, x=timescale_MinAge, color=TaxonomicClass2)) + 
  facet_grid(MinAge_mo~., scales="free") + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (generations)") + 
  geom_point(size=2, alpha=0.5) + scale_x_log10(breaks=seq(10,100,10)) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")

#results with shortened time series
ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=log10(timescale_MinAge), color=TaxonomicClass2, group=MainID)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("log10 Generation time (months)") + xlab("Time Series Length (log10 generations)") + 
  geom_point(data=gpdd_combo, size=1) + 
  geom_path(data=gpdd_combo, alpha=0.3) +
  geom_point(size=2) + 
  classic + labs(color="Taxonomic\nClass") + legalpha 
ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=log10(timescale_MinAge), color=mincisign, group=MainID)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("log10 Generation time (months)") + xlab("Time Series Length (log10 generations)") + 
  geom_point(data=gpdd_combo, size=1) + 
  geom_path(data=gpdd_combo, alpha=0.3, color="gray") +
  geom_point(size=2) + 
  classic + labs(color="Classification") + legalpha 
ggplot(gpdd_d, aes(y=E, x=log10(timescale_MinAge), color=log10(MinAge_mo), group=MainID)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("E") + xlab("Time Series Length (log10 generations)") + 
  geom_point(data=gpdd_combo, size=1) + 
  geom_path(data=gpdd_combo, alpha=0.3) +
  geom_point(size=2) + 
  classic + labs(color="log10\nGeneration\ntime (months)") + scale_color_viridis_c()
ggplot(gpdd_d, aes(y=minci_mo, x=log10(timescale_MinAge), color=E, group=MainID)) + 
  facet_grid(TaxonomicClass2~.) + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
  geom_hline(yintercept = 0) +
  geom_point(data=gpdd_combo, size=1, alpha=0.7) + 
  geom_path(data=gpdd_combo, alpha=0.3, color="gray") +
  geom_point(size=2, alpha=0.7) + 
  classic + labs(color="E") + scale_color_viridis_c()
ggplot(gpdd_d, aes(y=minci_mo, x=log10(timescale_MinAge), color=log10(MinAge_mo), group=MainID)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
  geom_hline(yintercept = 0) +
  geom_point(data=gpdd_combo, size=1, alpha=0.7) + 
  geom_path(data=gpdd_combo, alpha=0.3) +
  geom_point(size=2, alpha=0.7) + 
  classic + labs(color="log10\nGeneration\ntime (months)") + scale_color_viridis_c()
ggplot(gpdd_d, aes(y=minci_mo, x=log10(MinAge_mo), color=log10(timescale_MinAge), group=MainID)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("LE lower bound (per month)") + xlab("log10 Generation time (months)") + 
  geom_hline(yintercept = 0) +
  geom_point(data=gpdd_combo, size=1, alpha=0.7) + 
  geom_path(data=gpdd_combo, alpha=0.3, color="gray") +
  geom_point(size=2, alpha=0.7) + 
  classic + labs(color="Time Series Length\n(log10 generations)") + scale_color_viridis_c()

#select some individual series
gpdd_sub=arrange(gpdd_d, TaxonomicClass2, desc(datasetlength)) %>% 
  group_by(TaxonomicClass2)  %>% 
  slice(1:5)
gpdd_combo_sub=filter(gpdd_combo, MainID %in% unique(gpdd_sub$MainID))

#select only chaotic series
gpdd_sub=filter(gpdd_d, mincisign=="chaotic" & datasetlength>30 & !is.na(timescale_MinAge)) %>% arrange(desc(timescale_MinAge))
gpdd_combo_sub=filter(gpdd_combo, MainID %in% unique(gpdd_sub$MainID))

ggplot(gpdd_sub, aes(y=factor(MainID, levels = unique(gpdd_sub$MainID)), x=timescale_MinAge, fill=mincisign, group=MainID)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("MainID") + xlab("Time series length (generations)") + 
  geom_path(data=gpdd_combo_sub, size=2, aes(color=MinAge_mo)) +
  geom_point(data=gpdd_combo_sub, pch=21, size=2) + 
  geom_point(size=2, pch=21, aes(fill=mincisign)) + scale_x_log10(breaks=c(1,10,100,1000,10000)) +
  classic + scale_color_viridis_c(trans="log10")+ labs(fill="Classification", color="Generation\ntime (months)")  

ggplot(gpdd_sub, aes(y=E, x=log10(timescale_MinAge), color=log10(MinAge_mo), group=MainID)) + 
  facet_grid(TaxonomicClass2~.) + 
  ylab("E") + xlab("Time Series Length (log10 generations)") + 
  #geom_point(data=gpdd_combo_sub, size=1) + 
  geom_path(data=gpdd_combo_sub, alpha=0.7) +
  geom_point(size=2) + 
  classic + labs(color="log10\nGeneration\ntime (months)") + scale_color_viridis_c()
ggplot(gpdd_sub, aes(y=minci_mo, x=log10(timescale_MinAge), color=log10(MinAge_mo), group=MainID)) + 
  #facet_grid(TaxonomicClass2~.) + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
  geom_hline(yintercept = 0) +
  #geom_point(data=gpdd_combo_sub, size=1, alpha=0.7) + 
  geom_path(data=gpdd_combo_sub, alpha=0.7) +
  geom_point(size=2, alpha=0.7) + 
  classic + labs(color="log10\nGeneration\ntime (months)") + scale_color_viridis_c()

ggplot(gpdd_d, aes(y=minci_mo, x=log10(timescale_MinAge), color=log10(MinAge_mo), group=MainID)) + 
  facet_grid(TaxonomicClass2~.) + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
  geom_hline(yintercept = 0) +
  #geom_point(data=gpdd_combo, size=1, alpha=0.7) + 
  geom_path(data=gpdd_combo, alpha=0.7) +
  geom_point(size=2, alpha=0.7) + 
  classic + labs(color="log10\nGeneration\ntime (months)") + scale_color_viridis_c()
ggplot(gpdd_d, aes(y=minci_mo, x=datasetlength, color=log10(MinAge_mo), group=MainID)) + 
  facet_grid(TaxonomicClass2~.) + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (timesteps)") + 
  geom_hline(yintercept = 0) +
  geom_point(data=gpdd_combo, size=1, alpha=0.2) + 
  geom_path(data=gpdd_combo, alpha=0.2) +
  geom_point(size=2, alpha=0.7) + 
  classic + labs(color="log10\nGeneration\ntime (months)") + scale_color_viridis_c()

#linear models ####
m1=lm(minci_mo~E*log10(MinAge_mo), data=gpdd_d)
m1=lm(minci_mo~E*log10(timescale_MinAge), data=gpdd_d)
m1=lm(minci_mo~E*log10(timescale_MinAge)+TaxonomicClass2, data=filter(gpdd_d, TaxonomicClass2!="Other"))
m1=lm(minci_mo~E*log10(timescale_MinAge)*TaxonomicClass2, data=filter(gpdd_d, TaxonomicClass2!="Other"))
summary(m1)
anova(m1)



library(lme4)
library(car)
m1=lmer(minci_mo~E*log10(timescale_MinAge)+(1|TaxonomicClass2), data=droplevels(filter(gpdd_d, TaxonomicClass2!="Other")))
summary(m1)
Anova(m1)

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
