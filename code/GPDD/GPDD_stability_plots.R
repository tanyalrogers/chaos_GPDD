# Plots of GPDD chaos results
# Tanya Rogers

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

source("./code/Simulations/ggplot_themes_rogers.R")

#LE chaos detection results
gpdd_d=read.csv("./data/gpdd_ts_metadata.csv", stringsAsFactors = F)
gpdd_results=read.csv("./data/gpdd_results_smap.csv")
gpdd_d=left_join(gpdd_d, gpdd_results, by="MainID")
gpdd_combo=read.csv("./data/gpdd_results_truncation_smap.csv")

#results using other chaos detection methods
gpdd_other=read.csv("./data/gpdd_results_othermethods.csv")
gpdd_d=left_join(gpdd_d, gpdd_other, by="MainID")

gpdd_d$LEclass01=ifelse(gpdd_d$LEclass=="chaotic",1,0)
gpdd_d$RQA01=ifelse(gpdd_d$RQA=="chaotic",1,0)
gpdd_d$PE01=ifelse(gpdd_d$PE=="chaotic",1,0)

#general stats ####
#number of distinct taxa (138)
n_distinct(gpdd_d$TaxonID)
#LEs estimated (175)
length(which(!is.na(gpdd_d$LEmin)))
#prop pos LE
length(which(gpdd_d$LEmin>0.01))/length(which(!is.na(gpdd_d$LEmin)))
#prop pos LE 1d
length(which(gpdd_d$LEmin1d>0.01))/length(which(!is.na(gpdd_d$LEmin1d)))
#prop pos RQA
length(which(gpdd_d$RQA01==1))/length(which(!is.na(gpdd_d$RQA01)))
#prop pos PE
length(which(gpdd_d$PE01==1))/length(which(!is.na(gpdd_d$PE01)))

#classifications by taxon #####
signcountsnd=aggregate(LEmin~LEclass*TaxonomicClass3, data=gpdd_d, FUN=length) %>% 
  mutate(LEclass2=ifelse(LEclass=="not chaotic", LEmin*-1,LEmin), Econ="Free E")
signcounts1d=aggregate(LEmin1d~LEclass1d*TaxonomicClass3, data=gpdd_d, FUN=length) %>% 
  mutate(LEclass1d2=ifelse(LEclass1d=="not chaotic", LEmin1d*-1,LEmin1d), Econ="E = 1")
colnames(signcounts1d)<-colnames(signcountsnd)
signcounts=rbind(signcountsnd,signcounts1d); signcounts$Econ=factor(signcounts$Econ, levels=unique(signcounts$Econ))
ggplot(signcounts, aes(x=TaxonomicClass3, y=LEclass2)) + 
  facet_grid(.~Econ) + ylim(c(-60,20)) +
  geom_bar(aes(fill=LEclass), stat = "identity") +
  geom_hline(yintercept = 0)  + classic + xlabvert + ylab("Count") +
  labs(fill="Classification", x="Taxonomic Group") + removefacetbackground +
  scale_fill_brewer(palette = "Paired", direction = -1)
ggsave("./figures/classification.png", height = 3.5, width = 5)

#LE vs LE1d ####
ggplot(filter(gpdd_d, E>1), aes(y=LEmin1d, x=LEmin, fill=factor(E))) + 
  #facet_grid(predictable_ag~.) + 
  geom_point(size=3, pch=21, color="black", alpha=0.9) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + scale_fill_viridis_d() +
  classic + labs(x="LE, Free E", y="LE, E = 1", fill="Free E")
ggsave("./figures/LE1dvsnd.png", width = 4, height = 3)

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

#positive le vs mass #####
#comparison to Anderson and Gilooly 2020 (positive values only)
agle=read.csv("./data/AndersonGilloolyLEdata.csv", stringsAsFactors = F)
tglevels=c("Birds","Bony fishes", "Insects", "Mammals", "Phytoplankton","Zooplankton","Marine inverts", "Microbes")
agle$Tax4=factor(agle$TaxonomicClass3, levels = tglevels)
gpdd_d$Tax4=factor(gpdd_d$TaxonomicClass3, levels = tglevels)

lemass1=ggplot(filter(gpdd_d, LEmin_mo>0), aes(y=log10(LEmin_mo), x=log10(Mass_g))) + 
  ylab(expression(~log[10]~LE~(month^-1))) + xlab(expression(~log[10]~Mass~(g))) +
  geom_smooth(method="lm", se = F, color="black") +
  geom_smooth(data=agle, aes(y=log10(LE_mo)), method="lm", se = F, color="black", lty=2) +
  geom_point(aes(fill=Tax4, shape="GPDD"), color="black", alpha=0.9, stroke=0.5, size=3.5, show.legend = F) + 
  geom_point(data=agle, aes(y=log10(LE_mo), fill=Tax4, shape="AG2020"), stroke=0.5, size=2) +
  classic + labs(fill="Taxonomic\nGroup", shape="Source") + scale_fill_brewer(palette = "Dark2", drop=F) +
  scale_shape_manual(values = c(24,21)) + guides(fill = guide_legend(override.aes = list(size = 3.5, shape=21)),
                                                 shape = guide_legend(override.aes = list(size = c(2,3.5)))) +
  theme(axis.title = element_text(size=12), axis.text = element_text(size=10))
ggsave("./figures/mass_scaling.png", lemass1, width = 5.25, height = 4)

#classification and LE vs gen time, E ####
classgtE=ggplot(gpdd_d, aes(y=LEclass01, x=log10(MinAge_mo), color=E)) + 
  ylab("Proportion Chaotic") + xlab(expression(~log[10]~Generation~Time~(months))) + 
  geom_jitter(size=3, alpha=0.7, width=0, height=0.03) +
  stat_smooth(method="glm", method.args=list(family="binomial"), color="black") +
  classic + scale_color_viridis_c() 
LEmgtE=ggplot(gpdd_d, aes(y=LEmin_mo, x=log10(MinAge_mo), color=E)) + 
  ylab(expression(LE~(month^-1))) + xlab(expression(~log[10]~Generation~Time~(months))) + 
  geom_point(size=3, alpha=0.7) +
  geom_hline(yintercept = 0) +
  classic + scale_color_viridis_c() 
left=plot_grid(classgtE + theme(legend.position="none"), LEmgtE + theme(legend.position="none"), nrow = 1, labels=c("A","B"))
legend <- get_legend(LEmgtE + theme(legend.box.margin = margin(0, 0, 0, 3)))
plot_grid(left,legend, ncol = 2, rel_widths = c(1,0.1))
ggsave("./figures/gentime.png", width = 7, height = 3)

#classification by gen time by E (other methods)
ggplot(gpdd_d, aes(y=RQA01, x=log10(MinAge_mo), color=E)) + 
  ylab("Proportion Chaotic") + xlab(expression(~log[10]~Generation~Time~(months))) + 
  geom_point(size=3, alpha=0.7) +
  stat_smooth(method="glm", method.args=list(family="binomial")) +
  classic + scale_color_viridis_c() 
ggplot(gpdd_d, aes(y=PE01, x=log10(MinAge_mo), color=E)) + 
  ylab("Proportion Chaotic") + xlab(expression(~log[10]~Generation~Time~(months))) + 
  geom_point(size=3, alpha=0.7) +
  stat_smooth(method="glm", method.args=list(family="binomial")) +
  classic + scale_color_viridis_c() 

#LE vs R2 ####
r2class=ggplot(gpdd_d, aes(y=LEclass01, x=R2abund)) + 
  ylab("Proportion chaotic") + xlab(expression(R^2~'for'~abundance)) +
  geom_jitter(aes(fill=TaxonomicClass3), width=0, height=0.03, size=2.5, pch=21, alpha=0.9, color="black") +
  stat_smooth(method="glm", method.args=list(family="binomial"), color="black") +
  classic + labs(fill="Taxonomic\nGroup") +scale_fill_brewer(palette = "Dark2")
r2ab=ggplot(gpdd_d, aes(y=LEmin_mo, x=R2abund, fill=TaxonomicClass3)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab(expression(LE~(month^-1))) + xlab(expression(R^2~'for'~abundance)) +
  geom_hline(yintercept = 0) + #geom_vline(xintercept = 0.2, lty=2) + 
  geom_point(size=2.5, pch=21, alpha=0.9, color="black") +
  classic + labs(fill="Taxonomic\nGroup") + scale_fill_brewer(palette = "Dark2")
r2gr=ggplot(gpdd_d, aes(y=LEmin_mo, x=R2gr, fill=TaxonomicClass3)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab(expression(LE~(month^-1))) + xlab(expression(R^2~'for'~growth~rate)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0.2, lty=2) + 
  geom_point(size=2.5, pch=21, alpha=0.9, color="black") +
  classic + labs(fill="Taxonomic\nGroup") + scale_fill_brewer(palette = "Dark2")
r22=ggplot(gpdd_d, aes(y=R2abund, x=R2gr, fill=LEclass)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab(expression(R^2~'for'~abundance)) + xlab(expression(R^2~'for'~growth~rate)) +
  geom_point(size=2.5, pch=21, alpha=0.9, color="black") +
  geom_hline(yintercept = 0.2, lty=2) + geom_vline(xintercept = 0.2, lty=2) + 
  classic + labs(fill="Classification") + scale_fill_brewer(palette = "Paired", direction = -1)

r2plots=plot_grid(r2class + theme(legend.position="none"), r2ab + theme(legend.position="none"), nrow = 1, labels=c("A","B"))
legend <- get_legend(r2class + theme(legend.box.margin = margin(0, 0, 0, 3)))
plot_grid(r2plots,legend, ncol = 2, rel_widths = c(1,0.25))
ggsave("./figures/R2LE.png", width = 7, height = 3)

r2plots=plot_grid(r2ab + theme(legend.position="none"), r2gr + theme(legend.position="none"), align = 'vh',labels="AUTO")
legend <- get_legend(r2ab + theme(legend.box.margin = margin(0, 0, 0, 12)))
r2plots2=plot_grid(r22, legend, ncol=2, rel_widths = c(2.5,1), labels=c("C",""))
plot_grid(r2plots,r2plots2, nrow = 2)
ggsave("./figures/R2LE2.png", width = 6, height = 5)

#LE vs monotonic trend ####
mtclass=ggplot(gpdd_d, aes(y=LEclass01, x=monotonicR2)) + 
  ylab("Proportion chaotic") + xlab(expression(Monotonic~trend~R^2)) +
  geom_jitter(aes(fill=TaxonomicClass3), width=0, height=0.03, size=2.5, pch=21, alpha=0.9, color="black") +
  stat_smooth(method="glm", method.args=list(family="binomial"), color="black") +
  classic + labs(fill="Taxonomic\nGroup") +scale_fill_brewer(palette = "Dark2")
mtle=ggplot(gpdd_d, aes(y=LEmin_mo, x=monotonicR2, fill=TaxonomicClass3)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab(expression(LE~(month^-1))) + xlab(expression(Monotonic~trend~R^2)) +
  geom_point(size=2.5, pch=21, alpha=0.9) +
  geom_hline(yintercept = 0) +
  classic + labs(fill="Taxonomic\nClass") +scale_fill_brewer(palette = "Dark2")
mtplots=plot_grid(mtclass + theme(legend.position="none"), mtle + theme(legend.position="none"), nrow = 1, labels=c("A","B"))
legend <- get_legend(mtclass + theme(legend.box.margin = margin(0, 0, 0, 3)))
plot_grid(mtplots,legend, ncol = 2, rel_widths = c(1,0.25))
ggsave("./figures/monotonictrend.png", width = 7, height = 3)

#results with shortened chaotic series ####
gpdd_sub=filter(gpdd_d, LEclass=="chaotic" & datasetlength>30 & !is.na(timescale_MinAge)) %>% arrange(desc(timescale_MinAge))
gpdd_combo_sub=filter(gpdd_combo, MainID %in% unique(gpdd_sub$MainID)) 
ggplot(gpdd_sub, aes(y=factor(MainID, levels = unique(gpdd_sub$MainID)), x=timescale_MinAge, fill=LEclass, group=MainID)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab("MainID") + xlab("Generations Sampled") + 
  geom_path(data=gpdd_combo_sub, size=2, aes(color=MinAge_mo)) +
  geom_point(data=gpdd_combo_sub, pch=21, size=2) + 
  geom_point(size=2, pch=21, aes(fill=LEclass)) + scale_x_log10(breaks=c(1,10,100,1000,10000)) +
  classic + scale_color_viridis_c(trans="log10")+ scale_fill_manual(values = c("gray30","white")) +
  labs(fill="Classification", color="Generation\ntime (months)") + theme(axis.text.y = element_text(size=6)) 
ggsave("./figures/shortchaotic.png", width = 5, height = 5)

#change in classification full to 30
gpdd_combo_sub2=filter(gpdd_combo_sub, datasetlength==30) %>% left_join(select(gpdd_d, MainID, TaxonomicClass3, datasetlength), by="MainID")
gpdd_combo_sub2$LEcchange01=ifelse(gpdd_combo_sub2$LEclass=="chaotic",1,0)
gpdd_combo_sub2$changetsl=(gpdd_combo_sub2$datasetlength.y-gpdd_combo_sub2$datasetlength.x)/gpdd_combo_sub2$datasetlength.y
ggplot(gpdd_combo_sub2, aes(y=LEcchange01, x=log10(MinAge_mo), color=E)) + 
  ylab("Proportion remaining chaotic") + xlab(expression(~log[10]~Generation~Time~(months))) + 
  geom_jitter(size=3, alpha=0.7, width=0, height=0.03) +
  stat_smooth(method="glm", method.args=list(family="binomial"), color="black") +
  classic + scale_color_viridis_c() 
ggsave("./figures/shortchaotic_change.png", width = 4, height = 3)

mg=glm(LEcchange01~log10(MinAge_mo), data=gpdd_combo_sub2, family="binomial")
car::Anova(mg)
summary(mg)
