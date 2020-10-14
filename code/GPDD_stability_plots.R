# GPDD stability analysis - plots

library(ggplot2)
library(dplyr)
library(tidyr)
library(rEDM)
library(purrr)
library(cowplot)

source("./code/ggplot themes rogers.R")
#load(file = "./data/gpdd_results_update.Rdata")
gpdd_d=read.csv("./data/gpdd_ts_metadata.csv", stringsAsFactors = F)
gpdd_results=read.csv("./data/gpdd_results_smap.csv")
gpdd_d=left_join(gpdd_d, gpdd_results, by="MainID")
gpdd_combo=read.csv("./data/gpdd_results_truncation_smap.csv")
#write.csv(gpdd_d, "./data/gpdd_results_smap_steve.csv", row.names = F)

#general stats

#number of distinct taxa (138)
n_distinct(gpdd_d$TaxonID)
#gles estimated (175)
length(which(!is.na(gpdd_d$LEmin)))
#prop pos le
length(which(gpdd_d$LEmin>0.01))/length(which(!is.na(gpdd_d$LEmin)))
#prop pos le 1d
length(which(gpdd_d$LEmin1d>0.01))/length(which(!is.na(gpdd_d$LEmin1d)))
# old gle estimates
# #prop pos gles 
# length(which(gpdd_d$gle>0.01))/length(which(!is.na(gpdd_d$gle)))
# #prop pos gles 1d 
# length(which(gpdd_d$gle1d>0.01))/length(which(!is.na(gpdd_d$gle1d)))

#prop with any pos local le 
length(which(gpdd_d$LLE_proppos>0))/length(which(!is.na(gpdd_d$LLE_proppos)))
#prop with at least 0.5 pos local le 
length(which(gpdd_d$LLE_proppos>0.5))/length(which(!is.na(gpdd_d$LLE_proppos)))

#regression method
length(which(gpdd_d$LEregmin>0.01))/length(which(!is.na(gpdd_d$LEregmin)))
#sibly method
length(which(gpdd_d$LEsibly>0.01))/length(which(!is.na(gpdd_d$LEsibly)))

# #gles are correlated, although some outliers
# pairs(select(gpdd_d, gle1:gle5, gle))
# #timespan variable are correlated with each other, or their inverse
# pairs(select(gpdd_d, MinAge_mo, Lifespan_mo, timescale_MinAge, timescale_Lifespan, timestep_MinAge, timestep_Lifespan))

#gle is never higher than lle_avg
# plot(lle_avgbest~gle, data=gpdd_d, xlim=c(-1,1), ylim=c(-1,1))
# abline(a = 0, b = 1); abline(h=0); abline(v=0)
#gle vs LEmin (mostly correlated)
# plot(LEmin~gle, data=gpdd_d, xlim=c(-1,1), ylim=c(-1,1))
# abline(a = 0, b = 1); abline(h=0); abline(v=0)

#LE vs. LLE prop positive
# ggplot(gpdd_d, aes(x=LLE_proppos, fill=glesign)) + 
#   geom_histogram(boundary = 0, binwidth = 0.05) + classic
# ggplot(gpdd_d, aes(x=LLE_proppos, fill=LEclass)) + 
#   geom_histogram(boundary = 0, binwidth = 0.05) + classic

#reliability seems to be unrelated to predictability
# table(gpdd_d$Reliability, gpdd_d$predictable_ag)
# plot(gpdd_d$Reliability, gpdd_d$R2abund, ylim=c(-1,1))

#classifications by taxon (predictability) #####
signcountsnd=aggregate(LEmin~LEclass*TaxonomicClass3, data=gpdd_d, FUN=length) %>% 
  mutate(LEclass2=ifelse(LEclass=="not chaotic", LEmin*-1,LEmin), Econ="Free E") %>% 
  group_by(TaxonomicClass3) %>%  mutate(signprop=LEmin/sum(LEmin))
#for figure
signcountsnd=aggregate(LEmin~LEclass*TaxonomicClass3*predictable_ag, data=gpdd_d, FUN=length) %>% 
  mutate(LEclass2=ifelse(LEclass=="not chaotic", LEmin*-1,LEmin), Econ="Free E")
signcounts1d=aggregate(LEmin1d~LEclass1d*TaxonomicClass3*predictable_ag, data=gpdd_d, FUN=length) %>% 
  mutate(LEclass1d2=ifelse(LEclass1d=="not chaotic", LEmin1d*-1,LEmin1d), Econ="E = 1")
colnames(signcounts1d)<-colnames(signcountsnd)
signcounts=rbind(signcountsnd,signcounts1d); signcounts$Econ=factor(signcounts$Econ, levels=unique(signcounts$Econ))
ggplot(signcounts, aes(x=TaxonomicClass3, y=LEclass2)) + 
  #facet_grid(predictable_ag~., scales = "free_y") + #geom_jitter(alpha=0.5, size=3, height = 0, width = 0.1) +
  facet_grid(.~Econ) + ylim(c(-60,20)) +
  geom_bar(aes(fill=LEclass), stat = "identity") +
  geom_hline(yintercept = 0)  + classic + xlabvert + ylab("Count") +
  labs(fill="Classification", x="Taxonomic Group") + removefacetbackground +
  scale_fill_brewer(palette = "Paired", direction = -1)
ggsave("./figures/classification.png", height = 3.5, width = 5)

#classifications by E
# signcountsE=aggregate(LEmin~LEclass*E*TaxonomicClass3, data=gpdd_d, FUN=length) %>% 
#   mutate(LEclass2=ifelse(LEclass=="not chaotic", LEmin*-1,LEmin)) %>% 
#   group_by(E) %>% mutate(signprop=LEmin/sum(LEmin))
# ggplot(signcountsE, aes(x=factor(E), y=LEclass2)) + 
#   facet_grid(TaxonomicClass3~., scales = "free_y") +
#   geom_bar(aes(fill=LEclass), stat = "identity") +
#   geom_hline(yintercept = 0)  + classic + ylab("Count") +
#   labs(fill="Classification", x="E")
# ggplot(signcountsE, aes(x=factor(E), y=LEmin)) + 
#   facet_grid(TaxonomicClass3~., scales = "free_y") +
#   geom_bar(aes(fill=LEclass), stat = "identity") +
#   geom_hline(yintercept = 0)  + classic + ylab("Count") +
#   labs(fill="Classification", x="E")
# ggplot(signcountsE, aes(x=factor(E), y=signprop)) + 
#   geom_bar(aes(fill=LEclass), stat = "identity") +
#   classic + ylab("Proportion") + scale_y_continuous(expand = expand_scale(mult = c(0, 0))) +
#   labs(fill="Classification", x="E") 

#gle vs gle1d ####
ggplot(filter(gpdd_d, E>1), aes(y=LEmin1d, x=LEmin, fill=factor(E))) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=3, pch=21, color="black", alpha=0.9) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + scale_fill_viridis_d() +
  classic + labs(x="LE, Free E", y="LE, E = 1", fill="Free E")
ggsave("./figures/LE1dvsnd.png", width = 4, height = 3)

#histogram gle by taxon (predictability)
# ggplot(gpdd_d, aes(x=LEmin, fill=TaxonomicClass3)) + 
#   facet_grid(predictable_ag~.) + 
#   geom_histogram(boundary = 0.01, binwidth = 0.1) +
#   geom_vline(xintercept = 0) + theme_bw() + 
#   labs(x="LE lower bound", y="Count", fill="Taxonomic\nClass") 
# ggplot(gpdd_d, aes(x=LEmin1d, fill=TaxonomicClass3)) + 
#   #facet_grid(predictable_ag~.) + 
#   geom_histogram(boundary = 0, binwidth = 0.1) +
#   geom_vline(xintercept = 0) + theme_bw() 
# ggplot(gpdd_d, aes(x=LEmin_mo, fill=TaxonomicClass3)) + 
#   facet_grid(predictable_ag~.) + 
#   geom_histogram(boundary = 0, binwidth = 0.1) +
#   geom_vline(xintercept = 0) + theme_bw() 
# ggplot(gpdd_d, aes(x=LEmin_gen, fill=TaxonomicClass3)) + 
#   facet_grid(predictable_ag~.) + 
#   geom_histogram(boundary = 0, binwidth = 0.1) +
#   geom_vline(xintercept = 0) + theme_bw() + xlim(c(-2,2))
# #by sampling interval
# ggplot(gpdd_d, aes(x=LEmin, fill=TaxonomicClass3)) + 
#   facet_grid(SamplingInterval~.) + 
#   geom_histogram(boundary = 0, binwidth = 0.1) +
#   geom_vline(xintercept = 0) + theme_bw() 

#distribution of Es, taus ####
Es=ggplot(gpdd_d, aes(x=factor(E), fill=TaxonomicClass3)) + 
  #facet_grid(TaxonomicClass3~.) + 
  geom_bar(position = "stack") + xlab(expression(Embedding~dimension~(E))) +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Taxonomic\nClass") + scale_fill_brewer(palette = "Dark2")
taus=ggplot(gpdd_d, aes(x=factor(tau), fill=TaxonomicClass3)) + 
  #facet_grid(predictable_ag~.) + 
  geom_bar(position = "stack") + xlab(expression(Time~delay~(tau))) +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Taxonomic\nClass") + scale_fill_brewer(palette = "Dark2")

etau=plot_grid(Es + theme(legend.position="none"), taus + theme(legend.position="none"), align = 'vh',labels="AUTO")
legend <- get_legend(Es + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(etau, legend, rel_widths = c(3, 0.7))
ggsave("./figures/Etaudist.png", width = 8, height = 3.5)

#split by taxon
ggplot(gpdd_d, aes(x=factor(E), fill=TaxonomicClass3)) + 
  facet_grid(TaxonomicClass3~., scales="free_y") + 
  geom_bar(position = "stack") + xlab(expression(Embedding~dimension~(E))) +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Taxonomic\nClass") + scale_fill_brewer(palette = "Dark2")
ggplot(gpdd_d, aes(x=factor(tau), fill=factor(E))) + 
  #facet_grid(predictable_ag~.) + 
  geom_bar(position = "stack") + xlab("Best tau") +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Best E") 

#E vs time series length
# ggplot(gpdd_d, aes(y=E, x=datasetlength, color=TaxonomicClass3)) + 
#   #facet_grid(TaxonomicClass3~.) + 
#   ylab("E") + xlab("Time Series Length (timesteps)") + 
#   geom_point(size=2, alpha=0.4) + 
#   classic + labs(color="Taxonomic\nClass") + legalpha
# ggplot(gpdd_d, aes(y=E, x=timescale_MinAge*MinAge_mo/12, color=TaxonomicClass3)) + 
#   #facet_grid(TaxonomicClass3~.) + 
#   ylab("E") + xlab("Time Series Length (years)") + 
#   geom_point(size=2, alpha=0.4) + 
#   classic + labs(color="Taxonomic\nClass") + legalpha
# ggplot(gpdd_d, aes(y=E, x=log10(timescale_MinAge), color=TaxonomicClass3)) + 
#   #facet_grid(TaxonomicClass3~.) + 
#   ylab("E") + xlab("Time Series Length (log10 generations)") + 
#   geom_point(size=2, alpha=0.4) + 
#   classic + labs(color="Taxonomic\nClass") + legalpha
# #E vs gen time
# ggplot(gpdd_d, aes(y=E, x=log10(MinAge_mo), color=TaxonomicClass3)) + 
#   #facet_grid(TaxonomicClass3~.) + 
#   ylab("E") + xlab("log10 Generation Time (months)") + 
#   geom_point(size=2, alpha=0.4) + 
#   classic + labs(color="Taxonomic\nClass") + legalpha

#LE by E, theta
ggplot(gpdd_d, aes(y=LEmin, x=E)) + 
  xlab("E") + ylab("LE lower bound") +
  geom_point(aes(color=TaxonomicClass3), size=2, alpha=0.5) + 
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  geom_quantile(quantiles = c(0.20, 0.5, 0.80), formula=y~poly(x,2), size=2, alpha=0.5) +
  classic + labs(color="Taxonomic\nClass") + legalpha
ggplot(gpdd_d, aes(y=LEmin, x=factor(E), color=TaxonomicClass3)) + 
  #facet_grid(TaxonomicClass3~., scales="free_y") + 
  xlab("E") + ylab("LE lower bound") +
  geom_point(size=2, alpha=0.5) + 
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass") + legalpha
ggplot(gpdd_d, aes(y=LEmin, x=theta, color=TaxonomicClass3)) + 
  xlab("Theta") + ylab("LE lower bound") +
  geom_point(size=2, alpha=0.5) + 
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass") + legalpha

#age at maturity vs mass vs gens sampled ####
ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=log10(Mass_g), fill=log10(timescale_MinAge))) + 
  ylab(expression(log[10]~Age~at~Maturity~(months))) + 
  xlab(expression(log[10]~Mass~(g))) +
  geom_point(size=3, pch=21, alpha=0.9, color="black") +
  classic + labs(fill="log10\nGenerations\nsampled") + scale_fill_viridis_c()
ggsave("./figures/AgeMat_Mass_GenSampled.png", width = 4, height = 3)

# ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=log10(Lifespan_mo), fill=datasetlength)) + 
#   ylab("log10 Age at Maturity (months)") + xlab("log10 Lifespan (months)") + 
#   geom_point(size=2.5, pch=21, color="black") +
#   classic + labs(fill="Time Series\nLength")
# ggplot(gpdd_d, aes(y=log10(Lifespan_mo), x=log10(Mass_g), fill=datasetlength)) + 
#   ylab("log10 Lifespan (months)") + xlab("log10 Mass (g)") +
#   geom_point(size=2.5, pch=21, color="black") +
#   classic + labs(fill="Time Series\nLength")
# #gens/timespan vs gens/timestep
# ggplot(gpdd_d, aes(y=log10(timestep_MinAge), x=log10(MinAge_mo), fill=datasetlength)) + 
#   ylab("log10 Generations/Timestep") + xlab("log10 Age at Maturity (months)") +
#   geom_point(size=2.5, pch=21, color="black") +
#   classic + labs(fill="Time Series\nLength")
# ggplot(gpdd_d, aes(y=log10(timescale_MinAge), x=log10(timestep_MinAge), fill=datasetlength)) + 
#   ylab("log10 Generations sampled") + xlab("log10 Generations/Timestep") +
#   geom_point(size=2.5, pch=21, color="black") +
#   classic + labs(fill="Time Series\nLength")
# ggplot(gpdd_d, aes(y=log10(timescale_Lifespan), x=log10(timestep_Lifespan), fill=datasetlength)) + 
#   ylab("log10 Lifespans sampled") + xlab("log10 Lifespans/Timestep") +
#   geom_point(size=2.5, pch=21, color="black") +
#   classic + labs(fill="Time Series\nLength")

#gen time vs. ts length
# ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=datasetlength, color=TaxonomicClass3)) + 
#   #facet_grid(TaxonomicClass3~.) + 
#   ylab("log10 Generation time (months)") + xlab("Time Series Length (timesteps)") + 
#   geom_point(size=2, alpha=0.4) + 
#   classic + labs(color="Taxonomic\nClass") + legalpha
# ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=timescale_MinAge*MinAge_mo/12, color=TaxonomicClass3)) + 
#   #facet_grid(TaxonomicClass3~.) + 
#   ylab("log10 Generation time (months)") + xlab("Time Series Length (years)") + 
#   geom_point(size=2, alpha=0.4) + 
#   classic + labs(color="Taxonomic\nClass") + legalpha
# ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=log10(timescale_MinAge), color=TaxonomicClass3)) + 
#   #facet_grid(TaxonomicClass3~.) + 
#   ylab("log10 Generation time (months)") + xlab("Time Series Length (log10 generations)") + 
#   geom_point(size=2, alpha=0.4) + 
#   classic + labs(color="Taxonomic\nClass") + legalpha

#gle vs mass
ggplot(filter(gpdd_d, LEmin>0), aes(y=log10(LEmin_mo), x=log10(Mass_g), color=TaxonomicClass3)) +
  #facet_grid(predictable_ag~.) +
  #facet_grid(TaxonomicClass3~., scales="free_y") +
  ylab("LE lower bound (per month)") + xlab("log10 Mass (g)") +
  geom_point(size=2, alpha=0.9) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")
ggplot(filter(gpdd_d, LEmin>0), aes(y=log10(LEmin_mo), x=log10(MinAge_mo), color=TaxonomicClass3)) +
  #facet_grid(predictable_ag~.) +
  #facet_grid(TaxonomicClass3~., scales="free_y") +
  ylab("LE lower bound (per month)") + xlab("log10 Mass (g)") +
  geom_point(size=2, alpha=0.9) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")

#comparison to Anderson and Gilooly 2020 (positive values only)
agle=read.csv("./data/AndersonGilloolyLEdata.csv", stringsAsFactors = F)
tglevels=c("Birds","Bony fishes", "Insects", "Mammals", "Phytoplankton","Zooplankton","Marine inverts", "Microbes")
agle$Tax4=factor(agle$TaxonomicClass3, levels = tglevels)
gpdd_d$Tax4=factor(gpdd_d$TaxonomicClass3, levels = tglevels)

#positive le vs mass #####
lemass1=ggplot(filter(gpdd_d, LEmin>0), aes(y=log10(LEmin_mo), x=log10(Mass_g))) + 
  ylab(expression(~log[10]~LE~(month^-1))) + xlab(expression(~log[10]~Mass~(g))) +
  geom_smooth(method="lm", se = F, color="black") +
  geom_smooth(data=agle, aes(y=log10(LE_mo)), method="lm", se = F, color="black", lty=2) +
  geom_point(aes(fill=Tax4, shape="GPDD"), color="black", alpha=0.9, stroke=0.5, size=3.5, show.legend = F) + 
  geom_point(data=agle, aes(y=log10(LE_mo), fill=Tax4, shape="AG2020"), stroke=0.5, size=2) +
  classic + labs(fill="Taxonomic\nGroup", shape="Source") + scale_fill_brewer(palette = "Dark2", drop=F) +
  scale_shape_manual(values = c(24,21)) + guides(fill = guide_legend(override.aes = list(size = 3.5, shape=21)),
                                                 shape = guide_legend(override.aes = list(size = c(2,3.5))))
# lemass2=ggplot(filter(gpdd_d, LEmin>0), aes(y=log10(LEmin_gen), x=log10(Mass_g))) + 
#   ylab(expression(~log[10]~LE~(generation^-1))) + xlab(expression(~log[10]~Mass~(g))) +
#   geom_smooth(method="lm", se = F, color="black") +
#   geom_smooth(data=agle, aes(y=log10(LE*GT)), method="lm", se = F, color="black", lty=2, size=0.75) +
#   geom_hline(yintercept = log10(0.06), lty=2) +
#   geom_point(aes(fill=TaxonomicClass3, size="GPDD"), pch=21, color="black", alpha=0.9) + 
#   geom_point(data=agle, aes(y=log10(LE*GT), fill=factor(TaxonomicClass3), size="AG2020"), pch=21, stroke=0, show.legend = F) +
#   classic + labs(fill="Taxonomic\nClass", size="Source") + scale_fill_brewer(palette = "Dark2") +
#   scale_size_manual(values = c(2,3.5)) + guides(fill = guide_legend(override.aes = list(size = 3.5)))
# 
# lemass=plot_grid(lemass1 + theme(legend.position="none"), lemass2 + theme(legend.position="none"), align = 'vh',labels="AUTO")
# legend <- get_legend(lemass1 + theme(legend.box.margin = margin(0, 0, 0, 12)))
# plot_grid(lemass, legend, rel_widths = c(3, 0.7))
# ggsave("./figures/LEmass.png", width = 8, height = 3.5)

#slopes
summary(lm(log10(LEmin_mo)~log10(Mass_g), data=filter(gpdd_d, LEclass=="chaotic")))
summary(lm(log10(LEmin_mo)~log10(MinAge_mo), data=filter(gpdd_d, LEclass=="chaotic")))
summary(lm(log10(MinAge_mo)~log10(Mass_g), data=filter(gpdd_d, LEclass=="chaotic")))
summary(lm(log10(MinAge_mo)~log10(timescale_MinAge), data=filter(gpdd_d, LEclass=="chaotic")))
summary(lm(log10(LE_mo)~log10(Mass_g), data=agle))
summary(lm(log10(LEmin_gen)~log10(Mass_g), data=filter(gpdd_d, LEclass=="chaotic")))

#all LE vs mass (cropped)
ggplot(gpdd_d, aes(y=LEmin_gen, x=log10(Mass_g), color=TaxonomicClass3)) + 
  #facet_grid(predictable_ag~.) +
  #facet_grid(TaxonomicClass3~., scales="free_y") +
  ylab("LE lower bound (per generation)") + xlab("log10 Mass (g)") +
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass") + ylim(c(-.3,.3)) #ylim(c(-3,3))

#gle by gen time
ggplot(gpdd_d, aes(y=LEmin_mo, x=log10(MinAge_mo), color=TaxonomicClass3)) + 
  ylab("LE lower bound (per month)") + xlab("log10 Generation Time (months)") + 
  geom_point(size=2, alpha=0.4) +
  geom_hline(yintercept = 0) +
  classic + labs(color="Taxonomic\nClass") + legalpha
ggplot(gpdd_d, aes(y=LEmin_gen, x=log10(MinAge_mo), color=TaxonomicClass3)) + 
  ylab("LE lower bound (per generation)") + xlab("log10 Generation Time (months)") + 
  geom_point(size=2, alpha=0.4) +
  geom_hline(yintercept = 0) +
  classic + labs(color="Taxonomic\nClass") + ylim(c(-.3,.3)) + xlim(c(-2,0))
#gle by gen sampled
ggplot(gpdd_d, aes(y=LEmin_mo, x=log10(timescale_MinAge), color=TaxonomicClass3)) + 
  ylab("LE lower bound (per month)") + xlab("log10 Generations Sampled") + 
  geom_point(size=2, alpha=0.4) +
  geom_hline(yintercept = 0) +
  classic + labs(color="Taxonomic\nClass") + legalpha
ggplot(gpdd_d, aes(y=LEmin_gen, x=log10(timescale_MinAge), color=TaxonomicClass3)) + 
  ylab("LE lower bound (per generation)") + xlab("log10 Generations Sampled") + 
  geom_point(size=2, alpha=0.4) +
  geom_hline(yintercept = 0) +
  classic + labs(color="Taxonomic\nClass") + legalpha

#variance in gle by gen time, within species
ggplot(gpdd_d, aes(y=log10(LEvar_mo), x=log10(MinAge_mo), color=TaxonomicClass3)) + 
  ylab("log10 variance in LE (per month)") + xlab("log10 Generation Time (months)") + 
  geom_point(size=2, alpha=0.4) +
  classic + labs(color="Taxonomic\nClass") + legalpha
ggplot(gpdd_d, aes(y=log10(LEvar), x=log10(MinAge_mo), color=LEclass)) + 
  ylab("log10 variance in LE (per timestep)") + xlab("log10 Generation Time (months)") + 
  geom_point(size=2, alpha=0.4) +
  classic + labs(color="Classification") + legalpha

#variance gle (gen) by gen time
levar=gpdd_d %>% filter(!is.na(MinAge_mo)) %>% mutate(logGT=log10(MinAge_mo), GTbin=cut(logGT, breaks=c(seq(-2,1.5,0.5), max(logGT, na.rm=T)))) %>% 
  group_by(GTbin) %>% summarise(levar=var(LEmin_gen)) %>% mutate(GTmid=seq(-2,1.5,0.5)+0.25)
levarplot=ggplot(levar, aes(x=GTmid, y=log10(levar))) +
  geom_point(size=3, color="firebrick") +
  ylab(expression(~log[10]~var(LE)~(generation^-1))) + xlab(expression(~log[10]~Generation~Time~(months))) +
  classic + scale_x_continuous(limits = c(-2,2), breaks=c(seq(-2,2,1)))

#variance gle (gen) by gen sampled ####
levar=gpdd_d %>% filter(!is.na(timescale_MinAge)) %>% mutate(logGT=log10(timescale_MinAge), GTbin=cut(logGT, breaks=c(seq(0,3.5,0.5), max(logGT, na.rm=T)))) %>% 
  group_by(GTbin) %>% summarise(levar=var(LEmin_gen)) %>% mutate(GTmid=seq(0,3.5,0.5)+0.25)
levarplot=ggplot(levar, aes(x=GTmid, y=log10(levar))) +
  geom_point(size=3, color="firebrick") +
  ylab(expression(~log[10]~var(LE)~(generation^-1))) + xlab(expression(~log[10]~Generations~Sampled)) +
  classic + scale_x_continuous(limits = c(0.5,4), breaks=c(seq(0,4,1)))

#gle by gen time by E ####
LEmgtE=ggplot(gpdd_d, aes(y=LEmin_mo, x=log10(MinAge_mo), color=E)) + 
  ylab(expression(LE~(month^-1))) + xlab(expression(~log[10]~Generation~Time~(months))) + 
  geom_point(size=3, alpha=0.7) +
  geom_hline(yintercept = 0) +
  classic + scale_color_viridis_c() + 
  theme(legend.direction = "horizontal", legend.position = c(0,0),legend.justification = c(0,0), legend.background = element_rect(color="black"))
LEggtE=ggplot(gpdd_d, aes(y=LEmin_gen, x=log10(MinAge_mo), color=E)) + 
  ylab(expression(LE~(generation^-1))) + xlab(expression(~log[10]~Generation~Time~(months))) + 
  geom_point(size=3, alpha=0.7) +
  geom_hline(yintercept = 0) + 
  classic + scale_color_viridis_c() +
  theme(legend.direction = "horizontal", legend.position = c(0,0),legend.justification = c(0,0), legend.background = element_rect(color="black"))
# LEg2gtE=ggplot(gpdd_d, aes(y=log10(abs(LEmin_gen)), x=log10(MinAge_mo), color=LEclass)) +
#   ylab(expression(~log[10]~abs(LE)~(generation^-1))) + 
#   xlab(expression(~log[10]~Generation~Time~(months))) +  
#   geom_point(size=3, alpha=0.7) + labs(color="Classification") +
#   classic + scale_color_brewer(palette = "Paired", direction = -1)
# loglemass=ggplot(filter(gpdd_d, LEclass=="chaotic"), aes(y=log10(LEmin_mo), x=log10(Mass_g))) + 
#   ylab(expression(~log[10]~LE~(month^-1))) + xlab(expression(~log[10]~Mass~(g))) +
#   #geom_smooth(method="lm", se = F, color="black") +
#   geom_point(aes(fill=TaxonomicClass3), size=3, pch=21, color="black", alpha=0.7) + 
#   classic + labs(fill="Taxonomic\nClass") + scale_fill_brewer(palette = "Dark2")

left=plot_grid(LEmgtE, LEggtE, nrow = 2, labels=c("A","C"))
right=plot_grid(lemass1 + theme(legend.position="none"), levarplot, nrow = 2, labels=c("B","D"), axis = "btlr", align="hv")
legend <- get_legend(lemass1 + theme(legend.box.margin = margin(0, 0, 0, 3), legend.justification = "top"))
plot_grid(left,right,legend, ncol = 3, rel_widths = c(1,1,0.4))
ggsave("./figures/gentime.png", width = 8, height = 6)


# LEplots=plot_grid(LEmgtE + theme(legend.position="none"), LEg2gtE + theme(legend.position="none"), align = 'vh',labels=c("A","C"), ncol = 1)
# legend <- get_legend(LEmgtE + theme(legend.box.margin = margin(0, 0, 0, 12)))
# legend2 <- get_legend(LEg2gtE + theme(legend.box.margin = margin(0, 0, 0, 12)))
# legends <- plot_grid(legend2, legend, nrow = 1)
# LEplots2=plot_grid(LEggtE + theme(legend.position="none"), legends, ncol=1, labels=c("B",""))
# plot_grid(LEplots,LEplots2, nrow = 1)
# ggsave("./figures/gentime.png", width = 6, height = 5.1)
# 
# plot_grid(LEmgtE, LEggtE, LEg2gtE, nrow = 2, labels="AUTO", rel_widths = c(0.9,1))
# ggsave("./figures/gentime.png", width = 8, height = 6)

#gle by ts length
# ggplot(gpdd_d, aes(y=LEmin_mo, x=datasetlength, color=TaxonomicClass3)) + 
#   #facet_grid(TaxonomicClass3~., scales="free_y") + 
#   ylab("LE lower bound (per month)") + xlab("Time Series Length (timesteps)") +
#   geom_point(size=2, alpha=0.5) +
#   geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
#   classic + labs(color="Taxonomic\nClass")
# ggplot(gpdd_d, aes(y=LEmin_mo, x=log10(timescale_MinAge), color=TaxonomicClass3)) + 
#   #facet_grid(TaxonomicClass3~., scales="free_y") + 
#   ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
#   geom_point(size=2, alpha=0.5) +
#   geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
#   classic + labs(color="Taxonomic\nClass")
# #gle by ts length by E
# ggplot(gpdd_d, aes(y=LEmin_mo, x=log10(timescale_MinAge), color=E)) + 
#   #facet_grid(TaxonomicClass3~., scales="free_y") + 
#   ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
#   geom_point(size=2, alpha=0.5) +
#   geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
#   classic #+ labs(color="Taxonomic\nClass")
# #gle by ts length by gen time
# ggplot(gpdd_d, aes(y=LEmin_mo, x=log10(timescale_MinAge), color=log10(MinAge_mo))) + 
#   #facet_grid(TaxonomicClass3~., scales="free_y") + 
#   ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
#   geom_point(size=2, alpha=0.5) +
#   geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
#   classic + labs(color="log10\nGeneration\ntime (months)")

#R2 by time series length
# ggplot(gpdd_d, aes(y=R2abund, x=log10(timescale_MinAge), color=TaxonomicClass3)) + 
#   #facet_grid(predictable_ag~.) + 
#   ylab("R-squared") + xlab("Time Series Length (log10 Generations)") +
#   geom_point(size=2, alpha=0.5) +
#   geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
#   classic + labs(color="Taxonomic\nClass") + legalpha

#gle by R2 ####
r2ab=ggplot(gpdd_d, aes(y=LEmin_mo, x=R2abund, fill=TaxonomicClass3)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab(expression(LE~(month^-1))) + xlab(expression(R^2~'for'~abundance)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0.2, lty=2) + 
  geom_point(size=2.5, pch=21, alpha=0.9, color="black") +
  classic + labs(fill="Taxonomic\nClass") + scale_fill_brewer(palette = "Dark2")
r2gr=ggplot(gpdd_d, aes(y=LEmin_mo, x=R2gr, fill=TaxonomicClass3)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab(expression(LE~(month^-1))) + xlab(expression(R^2~'for'~growth~rate)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0.2, lty=2) + 
  geom_point(size=2.5, pch=21, alpha=0.9, color="black") +
  classic + labs(fill="Taxonomic\nClass") + scale_fill_brewer(palette = "Dark2")
r22=ggplot(gpdd_d, aes(y=R2abund, x=R2gr, fill=LEclass)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab(expression(R^2~'for'~abundance)) + xlab(expression(R^2~'for'~growth~rate)) +
  geom_point(size=2.5, pch=21, alpha=0.9, color="black") +
  geom_hline(yintercept = 0.2, lty=2) + geom_vline(xintercept = 0.2, lty=2) + 
  classic + labs(fill="Classification") + scale_fill_brewer(palette = "Paired", direction = -1)
r2plots=plot_grid(r2ab + theme(legend.position="none"), r2gr + theme(legend.position="none"), align = 'vh',labels="AUTO")
legend <- get_legend(r2ab + theme(legend.box.margin = margin(0, 0, 0, 12)))
r2plots2=plot_grid(r22, legend, ncol=2, rel_widths = c(2.5,1), labels=c("C",""))
plot_grid(r2plots,r2plots2, nrow = 2)
ggsave("./figures/R2LE.png", width = 6, height = 5)


#gle vs monotonic trend ####
ggplot(gpdd_d, aes(y=LEmin_mo, x=monotonicR2, fill=TaxonomicClass3)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab(expression(LE~(month^-1))) + xlab(expression(Monotonic~trend~R^2)) +
  geom_point(size=2.5, pch=21, alpha=0.9) +
  geom_hline(yintercept = 0) +
  classic + labs(fill="Taxonomic\nClass") +scale_fill_brewer(palette = "Dark2")
ggsave("./figures/monotonictrend.png", width = 4.5, height = 3)


#ts length vs monotonic trend
ggplot(gpdd_d, aes(y=monotonicR2, x=log10(timestep_MinAge), color=TaxonomicClass3)) + 
  #facet_grid(predictable_ag~.) + 
  ylab("Monotonic trend R-squared") + xlab("Time Series Length (log10 Generations)") +
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")
#R2 vs monotonic trend
ggplot(gpdd_d, aes(y=R2abund, x=monotonicR2, color=E)) + 
  facet_grid(TaxonomicClass3~.) + 
  ylab("R-squared") + xlab("Monotonic trend R-squared") +
  geom_point(size=2, alpha=0.5) +
  geom_abline(slope = 1, intercept = 0) + scale_color_viridis_c() +
  classic #+ labs(color="Taxonomic\nClass")
#E vs monotonic trend
ggplot(gpdd_d, aes(y=E, x=monotonicR2, color=LEclass)) + 
  facet_grid(TaxonomicClass3~.) + 
  ylab("E") + xlab("Monotonic trend R-squared") +
  geom_point(size=2, alpha=0.5) + #scale_color_viridis_c() +
  classic + labs(color="Classification")

#gle vs latitude
ggplot(gpdd_d, aes(y=LEmin, x=abs(LatDD), fill=TaxonomicClass3, color=LEclass)) + 
  facet_wrap(~TaxonomicClass3, nrow=2) +
  ylab("LE lower bound") + xlab("Latitude (abs. value)") + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Classification", fill="Taxonomic\nClass")
ggplot(gpdd_d, aes(y=LEmin_mo, x=abs(LatDD), fill=TaxonomicClass3, color=LEclass)) + 
  #facet_wrap(~TaxonomicClass3, nrow=2) +
  ylab("LE lower bound (per month)") + xlab("Latitude (abs. value)") + 
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Classification", fill="Taxonomic\nClass")
ggplot(gpdd_d, aes(y=LEmin_gen, x=abs(LatDD), fill=TaxonomicClass3, color=LEclass)) + 
  #facet_wrap(~TaxonomicClass3, nrow=2) +
  ylab("LE lower bound (per generation)") + xlab("Latitude (abs. value)") +
  geom_point(size=2, pch=21, stroke=1.5) +
  geom_hline(yintercept = 0) + scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Classification", fill="Taxonomic\nClass") +
  ylim(c(-3,3))

#gle vs growth rate
gpdd_d$meangr=map_dbl(gpdd_d$data_rescale, ~mean(.x$PopRescale_gr, na.rm=T))
gpdd_d$maxgr=map_dbl(gpdd_d$data_rescale, ~max(.x$PopRescale_gr, na.rm=T))
gpdd_d$maxgr_mo=gpdd_d$maxgr/timescale_mo(gpdd_d$SamplingInterval,1)

ggplot(gpdd_d, aes(y=LEmin_mo, x=meangr, color=TaxonomicClass3)) + 
  facet_grid(.~E) +
  ylab("LE lower bound (per month)") + xlab("Mean growth rate") + 
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")
ggplot(gpdd_d, aes(y=LEmin_mo, x=maxgr, color=TaxonomicClass3)) + 
  #facet_grid(.~E) +
  ylab("LE lower bound (per month)") + xlab("Max growth rate") + 
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")

ggplot(gpdd_d, aes(y=log10(maxgr_mo), x=log10(Mass_g), color=LEclass)) + 
  #facet_grid(.~E) +
  ylab("log10 Max growth rate (monthly)") + xlab("log10 Mass(g)") + 
  geom_point(size=2, alpha=0.5) +
  classic + labs(color="Classification")


#resolving tslength vs gen time
gpdd_d$log10_MinAge_mo_bin=cut(log10(gpdd_d$MinAge_mo), seq(from=-2, to=2.5, by=0.5))
ggplot(gpdd_d, aes(y=LEmin_mo, x=log10(timescale_MinAge), color=TaxonomicClass3)) + 
  facet_wrap(log10_MinAge_mo_bin~., nrow = 3, scales="free") + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")
ggplot(filter(gpdd_d, MinAge_mo %in% c(8,10,12,24)), aes(y=LEmin_mo, x=timescale_MinAge*MinAge_mo/12, color=TaxonomicClass3)) + 
  facet_grid(MinAge_mo~., scales="free") + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (years)") + 
  geom_point(size=2, alpha=0.5) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")
ggplot(filter(gpdd_d, MinAge_mo %in% c(8,10,12,24)), aes(y=LEmin_mo, x=timescale_MinAge, color=TaxonomicClass3)) + 
  facet_grid(MinAge_mo~., scales="free") + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (generations)") + 
  geom_point(size=2, alpha=0.5) + scale_x_log10(breaks=seq(10,100,10)) +
  geom_hline(yintercept = 0) + #scale_color_manual(values=c("red", "black", NA)) +
  classic + labs(color="Taxonomic\nClass")

#results with shortened time series
ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=log10(timescale_MinAge), color=TaxonomicClass3, group=MainID)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab("log10 Generation time (months)") + xlab("Time Series Length (log10 generations)") + 
  geom_point(data=gpdd_combo, size=1) + 
  geom_path(data=gpdd_combo, alpha=0.3) +
  geom_point(size=2) + 
  classic + labs(color="Taxonomic\nClass") + legalpha 
ggplot(gpdd_d, aes(y=log10(MinAge_mo), x=log10(timescale_MinAge), color=LEclass, group=MainID)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab("log10 Generation time (months)") + xlab("Time Series Length (log10 generations)") + 
  geom_point(data=gpdd_combo, size=1) + 
  geom_path(data=gpdd_combo, alpha=0.3, color="gray") +
  geom_point(size=2) + 
  classic + labs(color="Classification") + legalpha 
ggplot(gpdd_d, aes(y=E, x=log10(timescale_MinAge), color=log10(MinAge_mo), group=MainID)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab("E") + xlab("Time Series Length (log10 generations)") + 
  geom_point(data=gpdd_combo, size=1) + 
  geom_path(data=gpdd_combo, alpha=0.3) +
  geom_point(size=2) + 
  classic + labs(color="log10\nGeneration\ntime (months)") + scale_color_viridis_c()
ggplot(gpdd_d, aes(y=LEmin_mo, x=log10(timescale_MinAge), color=E, group=MainID)) + 
  facet_grid(TaxonomicClass3~.) + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
  geom_hline(yintercept = 0) +
  geom_point(data=gpdd_combo, size=1, alpha=0.7) + 
  geom_path(data=gpdd_combo, alpha=0.3, color="gray") +
  geom_point(size=2, alpha=0.7) + 
  classic + labs(color="E") + scale_color_viridis_c()
ggplot(gpdd_d, aes(y=LEmin_mo, x=log10(timescale_MinAge), color=log10(MinAge_mo), group=MainID)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
  geom_hline(yintercept = 0) +
  geom_point(data=gpdd_combo, size=1, alpha=0.7) + 
  geom_path(data=gpdd_combo, alpha=0.3) +
  geom_point(size=2, alpha=0.7) + 
  classic + labs(color="log10\nGeneration\ntime (months)") + scale_color_viridis_c()
ggplot(gpdd_d, aes(y=LEmin_mo, x=log10(MinAge_mo), color=log10(timescale_MinAge), group=MainID)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab("LE lower bound (per month)") + xlab("log10 Generation time (months)") + 
  geom_hline(yintercept = 0) +
  geom_point(data=gpdd_combo, size=1, alpha=0.7) + 
  geom_path(data=gpdd_combo, alpha=0.3, color="gray") +
  geom_point(size=2, alpha=0.7) + 
  classic + labs(color="Time Series Length\n(log10 generations)") + scale_color_viridis_c()

#select some individual series
gpdd_sub=arrange(gpdd_d, TaxonomicClass3, desc(datasetlength)) %>% 
  group_by(TaxonomicClass3)  %>% 
  slice(1:5)
gpdd_combo_sub=filter(gpdd_combo, MainID %in% unique(gpdd_sub$MainID))

#select only chaotic series
gpdd_sub=filter(gpdd_d, LEclass=="chaotic" & datasetlength>30 & !is.na(timescale_MinAge)) %>% arrange(desc(timescale_MinAge))
gpdd_combo_sub=filter(gpdd_combo, MainID %in% unique(gpdd_sub$MainID))

#shortened chaotic series ####
ggplot(gpdd_sub, aes(y=factor(MainID, levels = unique(gpdd_sub$MainID)), x=timescale_MinAge, fill=LEclass, group=MainID)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab("MainID") + xlab("Time series length (generations)") + 
  geom_path(data=gpdd_combo_sub, size=2, aes(color=MinAge_mo)) +
  geom_point(data=gpdd_combo_sub, pch=21, size=2) + 
  geom_point(size=2, pch=21, aes(fill=LEclass)) + scale_x_log10(breaks=c(1,10,100,1000,10000)) +
  classic + scale_color_viridis_c(trans="log10")+ scale_fill_manual(values = c("gray30","white")) +
  labs(fill="Classification", color="Generation\ntime (months)") + theme(axis.text.y = element_text(size=6)) 
ggsave("./figures/shortchaotic.png", width = 5, height = 5)

ggplot(gpdd_sub, aes(y=factor(MainID, levels = unique(gpdd_sub$MainID)), x=timescale_MinAge, fill=factor(E), group=MainID)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab("MainID") + xlab("Time series length (generations)") + 
  geom_path(data=gpdd_combo_sub, size=2, aes(color=MinAge_mo)) +
  geom_point(data=gpdd_combo_sub, pch=21, size=2) + 
  geom_point(size=2, pch=21, aes(fill=factor(E))) + scale_x_log10(breaks=c(1,10,100,1000,10000)) +
  classic + scale_color_viridis_c(trans="log10")+ scale_fill_viridis_d() +
  labs(fill="E", color="Generation\ntime (months)") + theme(axis.text.y = element_text(size=6)) 

#histogram of Es
ggplot(gpdd_sub, aes(x=factor(E), fill=TaxonomicClass3)) + 
  #facet_grid(TaxonomicClass3~.) + 
  geom_bar(position = "stack") + xlab(expression(Embedding~dimension~(E))) +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Taxonomic\nClass") + scale_fill_brewer(palette = "Dark2")
ggplot(filter(gpdd_combo_sub, datasetlength==30), aes(x=factor(E), fill=TaxonomicClass3)) + 
  #facet_grid(TaxonomicClass3~.) + 
  geom_bar(position = "stack") + xlab(expression(Embedding~dimension~(E))) +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Taxonomic\nClass") + scale_fill_brewer(palette = "Dark2")

#change in E and LE with short ts ####
shortsum=select(gpdd_sub, MainID, TaxonomicClass3, datasetlength, timescale_MinAge, MinAge_mo, Mass_g, E, LEmin_mo)
shortsum2=select(gpdd_combo_sub, MainID, datasetlength2=datasetlength, timescale_MinAge2=timescale_MinAge, E2=E, LEmin_mo2=LEmin_mo) %>% 
  filter(datasetlength2==30)
shortsum=left_join(shortsum, shortsum2)
shortsum=mutate(shortsum,
                dE=E2-E,
                dLEmin_mo=LEmin_mo2-LEmin_mo,
                dlogtimescale_MinAge=log10(timescale_MinAge2)-log10(timescale_MinAge)) 
d1=ggplot(shortsum) +
  geom_hline(yintercept = 0, lty=2) +
  labs(x="", y="Change in E") +
  geom_boxplot(aes(y=dE, x="")) + classic + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
d2=ggplot(shortsum) +
  geom_hline(yintercept = 0, lty=2) +
  labs(x="", y="Change in LE") +
  geom_boxplot(aes(y=dLEmin_mo, x="")) + classic + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
plot_grid(d1,d2,align = "hv", labels = "AUTO")
ggsave("./figures/changeE_LE.png", width = 4, height = 3)

ggplot(shortsum, aes(y=dE, x=dlogtimescale_MinAge, color=TaxonomicClass3)) +
  geom_hline(yintercept = 0, lty=2) +
  labs(x="Change in Generation Sampled", y="Change in E") +
  geom_point(size=3, alpha=0.5) + classic 
ggplot(shortsum, aes(y=LEmin_mo2, x=LEmin_mo, color=dE)) +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2) +
  labs(x="LE", y="LE short") +
  geom_point(size=3, alpha=0.5) + classic 
ggplot(shortsum, aes(y=dLEmin_mo, x=dE, color=TaxonomicClass3)) +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2) +
  labs(x="Change in E", y="Change in LE") +
  geom_point(size=3, alpha=0.5) + classic 


ggplot(gpdd_sub, aes(y=E, x=log10(timescale_MinAge), color=log10(MinAge_mo), group=MainID)) + 
  facet_grid(TaxonomicClass3~.) + 
  ylab("E") + xlab("Time Series Length (log10 generations)") + 
  #geom_point(data=gpdd_combo_sub, size=1) + 
  geom_path(data=gpdd_combo_sub, alpha=0.7) +
  geom_point(size=2) + 
  classic + labs(color="log10\nGeneration\ntime (months)") + scale_color_viridis_c()
ggplot(gpdd_sub, aes(y=LEmin_mo, x=log10(timescale_MinAge), color=log10(MinAge_mo), group=MainID)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
  geom_hline(yintercept = 0) +
  #geom_point(data=gpdd_combo_sub, size=1, alpha=0.7) + 
  geom_path(data=gpdd_combo_sub, alpha=0.7) +
  geom_point(size=2, alpha=0.7) + 
  classic + labs(color="log10\nGeneration\ntime (months)") + scale_color_viridis_c()

ggplot(gpdd_d, aes(y=LEmin_mo, x=log10(timescale_MinAge), color=log10(MinAge_mo), group=MainID)) + 
  facet_grid(TaxonomicClass3~.) + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (log10 generations)") + 
  geom_hline(yintercept = 0) +
  #geom_point(data=gpdd_combo, size=1, alpha=0.7) + 
  geom_path(data=gpdd_combo, alpha=0.7) +
  geom_point(size=2, alpha=0.7) + 
  classic + labs(color="log10\nGeneration\ntime (months)") + scale_color_viridis_c()
ggplot(gpdd_d, aes(y=LEmin_mo, x=datasetlength, color=log10(MinAge_mo), group=MainID)) + 
  facet_grid(TaxonomicClass3~.) + 
  ylab("LE lower bound (per month)") + xlab("Time Series Length (timesteps)") + 
  geom_hline(yintercept = 0) +
  geom_point(data=gpdd_combo, size=1, alpha=0.2) + 
  geom_path(data=gpdd_combo, alpha=0.2) +
  geom_point(size=2, alpha=0.7) + 
  classic + labs(color="log10\nGeneration\ntime (months)") + scale_color_viridis_c()

mean(shortsum$dLEmin_mo)
sd(shortsum$dLEmin_mo)/sqrt(length(shortsum$dLEmin_mo))

var(shortsum$LEmin_mo)
var(shortsum$LEmin_mo2)

#linear models ####
m1=lm(LEmin_mo~E*log10(MinAge_mo), data=gpdd_d)
m1=lm(LEmin_mo~E*log10(timescale_MinAge), data=gpdd_d)
m1=lm(LEmin_mo~E*log10(timescale_MinAge)+TaxonomicClass3, data=filter(gpdd_d, TaxonomicClass3!="Other"))
m1=lm(LEmin_mo~E*log10(timescale_MinAge)*TaxonomicClass3, data=filter(gpdd_d, TaxonomicClass3!="Other"))
summary(m1)
anova(m1)



library(lme4)
library(car)
m1=lmer(LEmin_mo~E*log10(timescale_MinAge)+(1|TaxonomicClass3), data=droplevels(filter(gpdd_d, TaxonomicClass3!="Other")))
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
