# Plots of GPDD chaos results

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

source("./code/Simulations/ggplot_themes.R")

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
#LEs estimated (172)
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
aggregate(LEmin~LEclass*TaxonomicClass3, data=gpdd_d, FUN=length) %>% 
  spread(LEclass,LEmin) %>% mutate(propchaotic=chaotic/(chaotic+`not chaotic`))
         
signcountsnd=aggregate(LEmin~LEclass*TaxonomicClass3, data=gpdd_d, FUN=length) %>% 
  mutate(LEclass2=ifelse(LEclass=="not chaotic", LEmin*-1,LEmin), Econ="Free E")
signcounts1d=aggregate(LEmin1d~LEclass1d*TaxonomicClass3, data=gpdd_d, FUN=length) %>% 
  mutate(LEclass1d2=ifelse(LEclass1d=="not chaotic", LEmin1d*-1,LEmin1d), Econ="E = 1")
colnames(signcounts1d)<-colnames(signcountsnd)
signcounts=rbind(signcountsnd,signcounts1d); signcounts$Econ=factor(signcounts$Econ, levels=unique(signcounts$Econ))
ggplot(signcounts, aes(x=TaxonomicClass3, y=LEclass2)) + 
  facet_grid(.~Econ) + ylim(c(-60,20)) +
  geom_bar(aes(fill=LEclass), stat = "identity", color="gray30") +
  geom_hline(yintercept = 0)  + classic + xlabvert + ylab("Count") +
  labs(fill="Classification", x="Taxonomic Group") + removefacetbackground +
  scale_fill_brewer(palette = "Paired", direction = -1)
ggsave("./figures/classification.png", height = 3.5, width = 5)
ggsave("./figures/Rogers_Fig2.pdf", height = 3.5, width = 5)

#LE vs LE1d ####
ggplot(filter(gpdd_d, E>1), aes(y=LEmean1d, x=LEmean, fill=factor(E))) + 
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
  labs(y="Count", fill="Taxonomic\nGroup") + scale_fill_brewer(palette = "Dark2")
taus=ggplot(gpdd_d, aes(x=factor(tau), fill=TaxonomicClass3)) + 
  geom_bar(position = "stack") + xlab(expression(Time~delay~(tau))) +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Taxonomic\nGroup") + scale_fill_brewer(palette = "Dark2")
etau=plot_grid(Es + theme(legend.position="none"), taus + theme(legend.position="none"), align = 'vh',labels="AUTO")
legend <- get_legend(Es + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(etau, legend, rel_widths = c(3, 0.7))
ggsave("./figures/Etaudist.png", width = 8, height = 3.5)

# E tau theta chaotic/notchaotic
Es=ggplot(gpdd_d, aes(x=factor(E), fill=TaxonomicClass3)) + 
  facet_grid(.~LEclass) + 
  geom_bar(position = "stack") + xlab(expression(Embedding~dimension~(E))) +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Taxonomic\nGroup") + scale_fill_brewer(palette = "Dark2") + 
  removefacetbackground
taus=ggplot(gpdd_d, aes(x=factor(tau), fill=TaxonomicClass3)) + 
  facet_grid(.~LEclass) + 
  geom_bar(position = "stack") + xlab(expression(Time~delay~(tau))) +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Taxonomic\nGroup") + scale_fill_brewer(palette = "Dark2") + 
  removefacetbackground
thetas=ggplot(gpdd_d, aes(x=factor(theta), fill=TaxonomicClass3)) + 
  facet_grid(.~LEclass) + 
  geom_bar(position = "stack") + xlab(expression(Nonlinearity~(theta))) +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Taxonomic\nGroup") + scale_fill_brewer(palette = "Dark2") + 
  removefacetbackground
plot_grid(Es,taus,thetas,nrow=3,labels = "AUTO")
ggsave("./figures/Etauthetadist_class.png", width = 8, height = 7)

#positive le vs mass #####

#comparison to Anderson and Gilooly 2020 (positive values only)
agle=read.csv("./data/AndersonGilloolyLEdata.csv", stringsAsFactors = F)
tglevels=c("Birds","Bony fishes", "Insects", "Mammals", "Phytoplankton","Zooplankton","Marine inverts", "Microbes")
agle$Tax4=factor(agle$TaxonomicClass3, levels = tglevels)
gpdd_d$Tax4=factor(gpdd_d$TaxonomicClass3, levels = tglevels)
agle$LEmin=NA

#supplemental lake time series
lakes=read.csv("./data/lakes_ts_metadata.csv",stringsAsFactors = F)
lakes2=read.csv("./data/lakes_results_smap.csv",stringsAsFactors = F)
lakes=lakes %>% filter(JLEsign=="chaotic") %>% left_join(lakes2)
lakes$Tax4=factor(lakes$TaxonomicClass3, levels = tglevels)

gpdd_pos=filter(gpdd_d, LEmin_mo>0) %>% select(Tax4,LEmin=LEmin_mo,LEmean=LEmean_mo, Mass_g) %>% mutate(Origin="field", Source="GPDD", stroke=0.5)
agle_pos=select(agle, Tax4, LEmin=LEmin, LEmean=LE_mo,Mass_g,Origin) %>% mutate(Source="AG2020", stroke=0.5)
lakes_pos=select(lakes, Tax4, LEmin=JLE, LEmean=JLEmean, Mass_g) %>% mutate(Origin="field",Source="Lakes", stroke=1.25)
posLE=rbind(gpdd_pos, agle_pos,lakes_pos) 
posLE2=rbind(gpdd_pos, agle_pos) 
#write.csv(posLE,"archive/posLE.csv",row.names = F)

lemass1=ggplot(posLE, aes(y=log10(LEmean), x=log10(Mass_g))) + 
  ylab(expression(~log[10]~LE~(month^-1))) + xlab(expression(~log[10]~Mass~(g))) +
  geom_smooth(data=posLE2, method="lm", se = F, color="black") +
  geom_point(aes(shape=Source,size=Source,fill=Tax4,stroke=stroke), color="black", alpha=0.75) + 
  classic + labs(fill="Taxonomic\nGroup", shape="Source") + scale_fill_brewer(palette = "Dark2", drop=F) +
  scale_color_brewer(palette = "Dark2", drop=F) +
  scale_shape_manual(values = c(24,21,22)) + scale_size_manual(values = c(2,3.5,3)) +
  guides(fill = guide_legend(override.aes = list(size = 3.5, shape=21))) +
  guides(shape = guide_legend(override.aes = list(stroke = c(0.5,0.5,1.25)))) +
  theme(axis.title = element_text(size=12), axis.text = element_text(size=10))
ggsave("./figures/mass_scaling.png", lemass1, width = 5.25, height = 4)
ggsave("./figures/Rogers_Fig4.pdf", lemass1, width = 5.25, height = 4)

lemass2=ggplot(posLE, aes(y=log10(LEmean), x=log10(Mass_g))) + 
  ylab(expression(~log[10]~LE~(month^-1))) + xlab(expression(~log[10]~Mass~(g))) +
  geom_smooth(data=posLE2, method="lm", se = F, color="black") +
  geom_linerange(aes(ymin=log10(LEmin),ymax=log10(LEmean), color=Tax4), alpha=1,show.legend = F) + 
  geom_point(aes(shape=Source,size=Source,fill=Tax4,stroke=stroke), color="black", alpha=0.75) + 
  classic + labs(fill="Taxonomic\nGroup", shape="Source") + scale_fill_brewer(palette = "Dark2", drop=F) +
  scale_color_brewer(palette = "Dark2", drop=F) +
  scale_shape_manual(values = c(24,21,22)) + scale_size_manual(values = c(2,3.5,3)-0.5) +
  guides(fill = guide_legend(override.aes = list(size = 3.5-0.5, shape=21))) +
  guides(shape = guide_legend(override.aes = list(stroke = c(0.5,0.5,1.25)))) +
  theme(axis.title = element_text(size=12), axis.text = element_text(size=10))
ggsave("./figures/mass_scaling_errorbar.png", lemass2, width = 5.25, height = 4)

#scaling exponent
summary(sexp<-lm(log10(LEmean)~log10(Mass_g),data=posLE2))

t1<-lm(log10(LEmean)~log10(Mass_g),data=gpdd_pos)
t2<-lm(log10(LEmean)~log10(Mass_g)+Tax4,data=gpdd_pos)
t3<-lm(log10(LEmean)~log10(Mass_g)*Tax4,data=gpdd_pos)
lmtest::lrtest(t1,t2)
lmtest::lrtest(t1,t3)
lmtest::lrtest(t2,t3)
anova(t3)

#lake predictions
lakepred=predict(sexp,newdata = lakes_pos)
ve=sum(residuals(sexp)^2)/sexp$df.residual
lakepredz=abs(lakepred-log10(lakes_pos$LEmean))/sqrt(ve)
lakepredz>1.96

#classification and LE vs gen time, E ####
classgtE=ggplot(gpdd_d, aes(y=LEclass01, x=log10(MinAge_mo), color=E)) + 
  ylab("Proportion Chaotic") + xlab(expression(~log[10]~Generation~Time~(months))) + 
  geom_jitter(size=3, alpha=0.7, width=0, height=0.03) +
  stat_smooth(method="glm", method.args=list(family="binomial"), color="black") +
  classic + scale_color_viridis_c() 
LEmgtE=ggplot(gpdd_d, aes(y=LEmean_mo, x=log10(MinAge_mo), color=E)) + 
  ylab(expression(LE~(month^-1))) + xlab(expression(~log[10]~Generation~Time~(months))) + 
  #geom_linerange(aes(ymin=LEmin_mo,ymax=LEmean_mo), alpha=1,show.legend = F) + 
  geom_point(size=3, alpha=0.7) +
  geom_hline(yintercept = 0) +
  classic + scale_color_viridis_c() 
left=plot_grid(classgtE + theme(legend.position="none"), LEmgtE + theme(legend.position="none"), nrow = 1, labels=c("A","B"))
legend <- get_legend(LEmgtE + theme(legend.box.margin = margin(0, 0, 0, 3)))
plot_grid(left,legend, ncol = 2, rel_widths = c(1,0.1))
ggsave("./figures/gentime.png", width = 7, height = 3)
ggsave("./figures/Rogers_Fig3.pdf", width = 7, height = 3)

mod=glm(LEclass01~log10(MinAge_mo), data=gpdd_d, family="binomial")
car::Anova(mod)

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

#E vs gen time
ggplot(gpdd_d, aes(y=E, x=log10(MinAge_mo), color=E)) + 
  ylab("Proportion Chaotic") + xlab(expression(~log[10]~Generation~Time~(months))) + 
  geom_point(size=3, alpha=0.7) +
  stat_smooth(method="lm") +
  classic + scale_color_viridis_c() 
cor(log10(gpdd_d$MinAge_mo),gpdd_d$E,use = "p")

#LE vs R2 ####
r2class=ggplot(gpdd_d, aes(y=LEclass01, x=R2abund)) + 
  ylab("Proportion chaotic") + xlab(expression(R^2~'for'~abundance)) +
  geom_jitter(aes(fill=TaxonomicClass3), width=0, height=0.03, size=2.5, pch=21, alpha=0.9, color="black") +
  stat_smooth(method="glm", method.args=list(family="binomial"), color="black") +
  classic + labs(fill="Taxonomic\nGroup") +scale_fill_brewer(palette = "Dark2")
r2ab=ggplot(gpdd_d, aes(y=LEmean_mo, x=R2abund, fill=TaxonomicClass3)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab(expression(LE~(month^-1))) + xlab(expression(R^2~'for'~abundance)) +
  geom_hline(yintercept = 0) + #geom_vline(xintercept = 0.2, lty=2) + 
  geom_point(size=2.5, pch=21, alpha=0.9, color="black") +
  classic + labs(fill="Taxonomic\nGroup") + scale_fill_brewer(palette = "Dark2")
r2gr=ggplot(gpdd_d, aes(y=LEmean_mo, x=R2gr, fill=TaxonomicClass3)) + 
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
r2plots1=plot_grid(r2plots,legend, ncol = 2, rel_widths = c(1,0.25))
#ggsave("./figures/R2LE.png", r2plots1, width = 7, height = 3)

r2plots=plot_grid(r2ab + theme(legend.position="none"), r2gr + theme(legend.position="none"), align = 'vh',labels="AUTO")
legend <- get_legend(r2ab + theme(legend.box.margin = margin(0, 0, 0, 12)))
r2plots2=plot_grid(r22, legend, ncol=2, rel_widths = c(2.5,1), labels=c("C",""))
plot_grid(r2plots,r2plots2, nrow = 2)
ggsave("./figures/R2LE2.png", width = 6, height = 5)

mod=glm(LEclass01~R2abund, data=gpdd_d, family="binomial")
car::Anova(mod)

#restrict to high R2 values
gpdd_sub=filter(gpdd_d,R2abund>=0.25)
length(which(gpdd_sub$LEmin>0.01))/length(which(!is.na(gpdd_sub$LEmin)))

#LE vs monotonic trend ####
mtclass=ggplot(gpdd_d, aes(y=LEclass01, x=monotonicR2)) + 
  ylab("Proportion chaotic") + xlab(expression(Monotonic~trend~R^2)) +
  geom_jitter(aes(fill=TaxonomicClass3), width=0, height=0.03, size=2.5, pch=21, alpha=0.9, color="black") +
  stat_smooth(method="glm", method.args=list(family="binomial"), color="black") +
  classic + labs(fill="Taxonomic\nGroup") +scale_fill_brewer(palette = "Dark2")
mtle=ggplot(gpdd_d, aes(y=LEmean_mo, x=monotonicR2, fill=TaxonomicClass3)) + 
  #facet_grid(TaxonomicClass3~.) + 
  ylab(expression(LE~(month^-1))) + xlab(expression(Monotonic~trend~R^2)) +
  geom_point(size=2.5, pch=21, alpha=0.9) +
  geom_hline(yintercept = 0) +
  classic + labs(fill="Taxonomic\nClass") +scale_fill_brewer(palette = "Dark2")
mtplots=plot_grid(mtclass + theme(legend.position="none"), mtle + theme(legend.position="none"), nrow = 1, labels=c("C","D"))
legend <- get_legend(mtclass + theme(legend.box.margin = margin(0, 0, 0, 3)))
mtplots1=plot_grid(mtplots,legend, ncol = 2, rel_widths = c(1,0.25))
#ggsave("./figures/monotonictrend.png", mtplots1, width = 7, height = 3)

#combined R2 and monotonic trend plot for paper
plot_grid(r2plots1,mtplots1,ncol=1)
ggsave("./figures/R2monotonicLE.png", width = 7, height = 5)
ggsave("./figures/Rogers_ED_Fig2.tiff",compression="lzw", width = 7, height = 5)

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

# histograms ####

#CV
breaks=pretty(gpdd_d$CV,10)
gpdd_d$CVclass=cut(gpdd_d$CV,breaks,labels = (breaks+diff(breaks[1:2])/2)[-length(breaks)])
signcountsCV=aggregate(LEmin~LEclass*CVclass, data=gpdd_d, FUN=length) %>% 
  mutate(LEclass2=ifelse(LEclass=="not chaotic", LEmin*-1,LEmin))
hCV=ggplot(signcountsCV, aes(x=CVclass, y=LEclass2)) + 
  geom_col(aes(fill=LEclass), width=1,color="gray30") +
  geom_hline(yintercept = 0)  + classic +  ylab("Count") +
  labs(fill="Classification", x="Coefficient of variation") + 
  scale_fill_brewer(palette = "Paired", direction = -1) +
  theme(legend.position = c(1,0), legend.justification = c(1,0), legend.background = element_blank())

#R2
breaks=pretty(gpdd_d$R2abund,10)
gpdd_d$R2class=cut(gpdd_d$R2abund,breaks,labels = (breaks+diff(breaks[1:2])/2)[-length(breaks)])
signcountsR2=aggregate(LEmin~LEclass*R2class, data=gpdd_d, FUN=length) %>% 
  mutate(LEclass2=ifelse(LEclass=="not chaotic", LEmin*-1,LEmin))
hR2=ggplot(signcountsR2, aes(x=R2class, y=LEclass2)) + 
  geom_col(aes(fill=LEclass), width=1,color="gray30", show.legend = F) +
  geom_hline(yintercept = 0)  + classic +  ylab("Count") +
  labs(fill="Classification", x=expression(R^2~'for'~abundance)) + 
  scale_fill_brewer(palette = "Paired", direction = -1) 

#theta
signcountstheta=aggregate(LEmin~LEclass*theta, data=gpdd_d, FUN=length) %>% 
  mutate(LEclass2=ifelse(LEclass=="not chaotic", LEmin*-1,LEmin))
htheta=ggplot(signcountstheta, aes(x=factor(theta), y=LEclass2)) + 
  ylim(c(-61,20)) +
  geom_col(aes(fill=LEclass), color="gray30", show.legend = F) +
  geom_hline(yintercept = 0)  + classic +  ylab("Count") +
  labs(fill="Classification", x=expression(Nonlinearity~(theta))) + removefacetbackground +
  scale_fill_brewer(palette = "Paired", direction = -1) 

#monotonic trend
breaks=pretty(gpdd_d$monotonicR2,10)
gpdd_d$monoclass=cut(gpdd_d$monotonicR2,breaks,labels = (breaks+diff(breaks[1:2])/2)[-length(breaks)])
signcountsmono=aggregate(LEmin~LEclass*monoclass, data=gpdd_d, FUN=length) %>% 
  mutate(LEclass2=ifelse(LEclass=="not chaotic", LEmin*-1,LEmin))
hmono=ggplot(signcountsmono, aes(x=monoclass, y=LEclass2)) + 
  geom_col(aes(fill=LEclass), width=1,color="gray30", show.legend = F) +
  geom_hline(yintercept = 0)  + classic +  ylab("Count") +
  labs(fill="Classification", x=expression(Monotonic~trend~R^2)) + 
  scale_fill_brewer(palette = "Paired", direction = -1)

plot_grid(hCV,hR2,htheta,hmono,ncol = 2, align = "vh", labels = "AUTO")
ggsave("./figures/histograms.png", width = 8, height = 6)
ggsave("./figures/Rogers_Fig1.pdf", width = 8, height = 6)

# restrict ts length and aggregate by site ####
# for site level, classify using median and max

gpdd_sub50=filter(gpdd_d,ndatapoints>=50)
gpdd_sub70=filter(gpdd_d,ndatapoints>=70)
gpdd_sub70b=filter(gpdd_d,ndatapoints>=70 & timestep_MinAge<4)

tstable=data.frame(subset=c("all","50","70","70b"),nseries=NA,nsites=NA,pindiv=NA,psitemax=NA,psitemed=NA)

tstable$nseries[1]=length(which(!is.na(gpdd_d$LEmin)))
tstable$nseries[2]=length(which(!is.na(gpdd_sub50$LEmin)))
tstable$nseries[3]=length(which(!is.na(gpdd_sub70$LEmin)))
tstable$nseries[4]=length(which(!is.na(gpdd_sub70b$LEmin)))

tstable$pindiv[1]=length(which(gpdd_d$LEmin>0.01))/length(which(!is.na(gpdd_d$LEmin)))
tstable$pindiv[2]=length(which(gpdd_sub50$LEmin>0.01))/length(which(!is.na(gpdd_sub50$LEmin)))
tstable$pindiv[3]=length(which(gpdd_sub70$LEmin>0.01))/length(which(!is.na(gpdd_sub70$LEmin)))
tstable$pindiv[4]=length(which(gpdd_sub70b$LEmin>0.01))/length(which(!is.na(gpdd_sub70b$LEmin)))

tstable$nsites[1]=n_distinct(gpdd_d$ExactName)
tstable$nsites[2]=n_distinct(gpdd_sub50$ExactName)
tstable$nsites[3]=n_distinct(gpdd_sub70$ExactName)
tstable$nsites[4]=n_distinct(gpdd_sub70b$ExactName)

table(gpdd_d$ExactName)
gpdd_avg=gpdd_d %>% group_by(ExactName) %>% 
  summarise(LEmedian=median(LEmin),
            LEmax=max(LEmin))
gpdd_avg50=gpdd_sub50 %>% group_by(ExactName) %>% 
  summarise(LEmedian=median(LEmin),
            LEmax=max(LEmin))
gpdd_avg70=gpdd_sub70 %>% group_by(ExactName) %>% 
  summarise(LEmedian=median(LEmin),
            LEmax=max(LEmin))
gpdd_avg70b=gpdd_sub70b %>% group_by(ExactName) %>% 
  summarise(LEmedian=median(LEmin),
            LEmax=max(LEmin))

tstable$psitemed[1]=length(which(gpdd_avg$LEmedian>0.01))/length(which(!is.na(gpdd_avg$LEmedian)))
tstable$psitemax[1]=length(which(gpdd_avg$LEmax>0.01))/length(which(!is.na(gpdd_avg$LEmax)))
tstable$psitemed[2]=length(which(gpdd_avg50$LEmedian>0.01))/length(which(!is.na(gpdd_avg50$LEmedian)))
tstable$psitemax[2]=length(which(gpdd_avg50$LEmax>0.01))/length(which(!is.na(gpdd_avg50$LEmax)))
tstable$psitemed[3]=length(which(gpdd_avg70$LEmedian>0.01))/length(which(!is.na(gpdd_avg70$LEmedian)))
tstable$psitemax[3]=length(which(gpdd_avg70$LEmax>0.01))/length(which(!is.na(gpdd_avg70$LEmax)))
tstable$psitemed[4]=length(which(gpdd_avg70b$LEmedian>0.01))/length(which(!is.na(gpdd_avg70b$LEmedian)))
tstable$psitemax[4]=length(which(gpdd_avg70b$LEmax>0.01))/length(which(!is.na(gpdd_avg70b$LEmax)))

tstable[,4:6]=round(tstable[,4:6],2)
tstable

# corrected site level prob chaos #### 
gpdd_sitetax=gpdd_d %>% group_by(ExactName) %>% 
  summarise(sitetax=paste0(sort(unique(TaxonomicClass3)),collapse = ", "))
gpdd_prop=gpdd_d %>% group_by(ExactName) %>% 
  summarise(nseries=n(),
            propchaotic=length(which(LEclass=="chaotic"))/length(LEclass),
            propcorrected=min(max(0,((propchaotic-0.04)/(0.71-0.04))),1)) %>% 
  left_join(gpdd_sitetax) %>% 
  mutate(sitetax=factor(sitetax,levels=c("Birds","Bony fishes","Insects","Mammals",
                                         "Phytoplankton","Phytoplankton, Zooplankton",
                                         "Birds, Mammals")))

length(which(gpdd_prop$propcorrected>0.5))/length(gpdd_prop$propcorrected)
length(which(gpdd_prop$propcorrected>0.8))/length(gpdd_prop$propcorrected)

p1=ggplot(gpdd_prop, aes(x=nseries, fill=sitetax)) + 
  geom_histogram(position = "stack",binwidth = 1,boundary=0.5) + xlab("Number of series per location") +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(y="Count", fill="Taxonomic Group") + scale_fill_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = seq(0,35,5), expand = expand_scale(mult = c(0, 0.05)))

p2=ggplot(gpdd_prop, aes(x=propcorrected, fill=sitetax)) + 
  geom_histogram(position = "stack",binwidth = 0.2,boundary=0) + xlab("Probability that 'system' is chaotic") +
  classic + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  labs(y="Count", fill="Taxonomic Group") + scale_fill_brewer(palette = "Dark2")
#ggsave("./figures/siteprobchaos.png", width = 5.5, height = 3.5)

pcplots=plot_grid(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), nrow = 1, labels=c("A","B"))
legend <- get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 3), legend.text = element_text(size=8)))
plot_grid(pcplots,legend, ncol = 2, rel_widths = c(1,0.3))
ggsave("./figures/histo_siteprobchaos.png", width = 8, height = 3)
ggsave("./figures/Rogers_ED_Fig3.tiff", compression="lzw", width = 8, height = 3)

# prob chaos by taxon
gpdd_avg = gpdd_avg %>% left_join(gpdd_sitetax)
gpdd_avg %>% group_by(sitetax) %>% 
  summarise(n=n(), chaotic = length(which(LEmedian>0.01)), propchaotic=chaotic/n) %>% 
  filter(n>3) %>% as.data.frame()

#lakes
lakes2 %>% group_by(Site) %>% 
  summarise(LEmedian=median(JLE),
            nseries=n(),
            propchaotic=length(which(JLEsign=="chaotic"))/nseries,
            propcorrected=min(max(0,((propchaotic-0.04)/(0.71-0.04))),1))

# histograms for sites with >8 time series ####

lehisto=filter(gpdd_d, ExactName %in% c("Port Erin Bay","William Trelease Woods","Hinkley Point"))
lehisto$TaxonomicClass3=factor(lehisto$TaxonomicClass3,levels = sort(unique(gpdd_d$TaxonomicClass3)))
ggplot(lehisto,aes(x=LEmean_mo,fill=TaxonomicClass3)) +
  facet_wrap(.~ExactName, scale="free") +
  geom_histogram(bins=25) +
  geom_vline(aes(xintercept=0)) +
  classic + labs(x=expression(LE~(month^-1)),y="Count", fill="Taxonomic\nGroup") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  scale_fill_brewer(palette = "Dark2", drop=F) +
  removefacetbackground
ggsave("./figures/lehistos.png", width = 8, height = 3)
ggsave("./figures/Rogers_ED_Fig4.tiff", compression="lzw", width = 8, height = 3)

# sample time series plots ####

gpdd_d_ts=read.csv("./data/gpdd_timeseries.csv",stringsAsFactors = F)

set.seed(111)
gpdd_samp=select(gpdd_d,MainID,TaxonomicClass3,LEclass,TaxonName,CommonName) %>% 
  group_by(TaxonomicClass3,LEclass) %>% 
  mutate(keep=MainID %in% sample(MainID,size = 1,replace = F)) %>% 
  filter(keep==T) %>% ungroup() %>% 
  # mutate(label=paste0(CommonName," (", TaxonName, ") - ",LEclass)) %>% 
  # mutate(label=paste(MainID, CommonName,"-",LEclass)) %>% 
  mutate(label=paste0(CommonName, " (",MainID, ")")) %>% 
  arrange(TaxonomicClass3,LEclass) %>% 
  mutate(label=factor(label,levels=unique(label)))

gpdd_d_ts_sub=right_join(gpdd_d_ts,gpdd_samp) 

ggplot(gpdd_d_ts_sub, aes(x=SeriesStep, y=PopRescale_log)) +
  facet_wrap(.~label, scales="free", ncol=2) +
  geom_line() + geom_point() +
  classic + removefacetbackground +
  labs(x="Timestep",y="log Scaled Abundance")
ggsave("figures/sample_timeseries.png",width = 6, height = 7.5)
ggsave("figures/Rogers_ED_Fig5.tiff", compression="lzw",width = 6, height = 7.5)
