# Plot simuation results

library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")

source("./code/ggplot themes rogers.R")

#load results
load("./data/sims_results_update.Rdata")

#reclass noise level for stochastic ts
sims_d$NoiseLevel2=ifelse(sims_d$NoiseLevel==0, 0.01, sims_d$NoiseLevel)

#class LEs
sims_d$gle_class=ifelse(sims_d$gle>0.01, "chaotic", "not chaotic")
sims_d$LEshift_class=ifelse(sims_d$minci>0.01, "chaotic", "not chaotic")
sims_d$LEshift_class.05=ifelse(sims_d$minci>0.05, "chaotic", "not chaotic")

sims_summary=sims_d %>% select(-data) %>% 
  group_by(Classification,NoiseLevel2,TSlength, Model) %>% 
  summarize(gle_pp=length(which(gle>0.0))/length(gle),
            gle_pp.01=length(which(gle>0.01))/length(gle),
            gle_pp.05=length(which(gle>0.05))/length(gle),
            LEshift_pp=length(which(minci>0))/length(minci),
            LEshift_pp.01=length(which(minci>0.01))/length(minci),
            LEshift_pp.05=length(which(minci>0.05))/length(minci))

#plots

#proportions, individual models
pdf("SimClassLE.pdf", width = 12, height = 7)
ggplot(sims_summary, aes(x=Model, y=LEshift_pp.01, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel2) + 
  geom_bar(stat = "identity", color="black") + classic + xlabvert + ylab("Proportion Classified Chaotic") +
  labs(fill="True\nModel\nDynamics") + 
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
dev.off()
#different arrangement
ggplot(sims_summary, aes(x=NoiseLevel2, y=LEshift_pp.01, fill=Classification)) +
  facet_grid(TSlength~Classification) + 
  geom_bar(stat = "identity", position = position_dodge2(), show.legend = F) + theme_bw() + xlabvert


#plot values

#all models together
ggplot(sims_d, aes(x=NoiseLevel2, y=gle, color=Classification)) +
  facet_grid(TSlength~., scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.05), alpha=0.4) + 
  theme_bw() + xlab("Noise Level") + legalpha
#violin plots
ggplot(sims_d, aes(x=factor(NoiseLevel2), y=gle, fill=Classification)) +
  facet_grid(TSlength~., scales = "free_y") + geom_hline(yintercept = 0) +
  geom_violin(position = position_dodge(0.9), scale = "width") + 
  theme_bw() + xlab("Noise Level")

#individual models
ggplot(sims_d, aes(x=Model, y=gle, color=Classification)) +
  facet_grid(TSlength~NoiseLevel2, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02), alpha=0.3) + theme_bw() + xlabvert
ggplot(sims_d, aes(x=Model, y=minci, color=Classification)) +
  facet_grid(TSlength~NoiseLevel2, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02), alpha=0.3) + theme_bw() + xlabvert

#R squared, E, theta
ggplot(sims_d, aes(x=Model, y=R2best, color=Classification)) +
  facet_grid(TSlength~NoiseLevel2, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02), alpha=0.3) + theme_bw() + xlabvert + legalpha + ylab("R-squared")
ggplot(sims_d, aes(x=Model, y=Ebest, color=Classification)) +
  facet_grid(TSlength~NoiseLevel2) + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02), alpha=0.3) + theme_bw() + xlabvert + legalpha
ggplot(sims_d, aes(x=Model, y=taubest, color=Classification)) +
  facet_grid(TSlength~NoiseLevel2) + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02), alpha=0.3) + theme_bw() + xlabvert + legalpha
ggplot(sims_d, aes(x=Model, y=thetabest, color=Classification)) +
  facet_grid(TSlength~NoiseLevel2) + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02), alpha=0.3) + theme_bw() + xlabvert + legalpha
#E vs tau
ggplot(sims_d, aes(x=Ebest, y=taubest, size=TSlength, color=NoiseLevel2)) +
  facet_wrap(Model~., nrow = 3, scales = "free_y") + geom_hline(yintercept = 1) +
  geom_point(alpha=0.5) + #geom_boxplot() + 
  theme_bw() 

#histograms of LE distributions all model
tslengths=unique(sims_d$TSlength)
for(i in 1:length(tslengths)) {
  print(ggplot(filter(sims_d, TSlength==tslengths[i]), aes(x=minci, fill=Classification)) +
          facet_grid(Classification~NoiseLevel2, scales = "free_y") + geom_hline(yintercept = 0) +
          geom_histogram(boundary = 0.01, binwidth = 0.1, show.legend = F) + geom_vline(xintercept = 0) +
          theme_bw() + ggtitle(paste("TSlength =", tslengths[i])) + xlim(c(-1,1)))
}
#freq polygons (not as good)
ggplot(sims_d, aes(x=minci, color=Classification)) +
  facet_grid(TSlength~NoiseLevel2, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_freqpoly(boundary = 0.01, binwidth = 0.1, size=1) + geom_vline(xintercept = 0) +
  theme_bw() + xlim(c(-1,1))


#proportion correct classifications
#LEshift 0.01
sims_summary2=sims_d %>% select(-data) %>% 
  group_by(Classification,LEshift_class,NoiseLevel2,TSlength) %>% 
  summarize(n=n()) %>% ungroup() %>% 
  complete(Classification, nesting(LEshift_class,NoiseLevel2,TSlength), 
           fill=list(n=0)) %>% 
  group_by(NoiseLevel2,TSlength, Classification) %>% 
  mutate(proportion=n/sum(n))

ggplot(sims_summary2, aes(x=LEshift_class, y=Classification, fill=proportion)) +
  facet_grid(TSlength~NoiseLevel2) + 
  geom_tile(stat = "identity") + classic +
  geom_text(aes(label=round(proportion,2)), color="white") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  ylab("True Dynamics") + xlab("Classification") +
  labs(fill="Proportion") 

#LEshift 0.05
sims_summary3=sims_d %>% select(-data) %>% 
  group_by(Classification,LEshift_class.05,NoiseLevel2,TSlength) %>% 
  summarize(n=n()) %>% ungroup() %>% 
  complete(Classification, nesting(LEshift_class.05,NoiseLevel2,TSlength), 
           fill=list(n=0)) %>% 
  group_by(NoiseLevel2,TSlength, Classification) %>% 
  mutate(proportion=n/sum(n))

ggplot(sims_summary3, aes(x=LEshift_class.05, y=Classification, fill=proportion)) +
  facet_grid(TSlength~NoiseLevel2) + 
  geom_tile(stat = "identity") + classic + 
  geom_text(aes(label=round(proportion,2)), color="white") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  ylab("True Dynamics") + xlab("Classification") +
  labs(fill="Proportion")

#regression method
ggplot(sims_d, aes(x=NoiseLevel, y=LEreg, color=Classification)) +
  facet_grid(TSlength~.) + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02)) + theme_bw()

ggplot(sims_d, aes(x=Model, y=LEreg, color=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point() + #geom_boxplot() + 
  theme_bw() + xlabvert

