# Plot simuation results

library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")

source("./code/ggplot themes rogers.R")

#test data ####

#load results
load("./data/sims_results_update.Rdata")

#reclass noise level for stochastic ts
sims_d$NoiseLevel2=ifelse(sims_d$NoiseLevel==0, 0.01, sims_d$NoiseLevel)

#class LEs
sims_d$LEsign=ifelse(sims_d$minci>0.01, "chaotic", "not chaotic")

sims_summary=sims_d %>% select(-data) %>% 
  group_by(Classification,NoiseLevel2,TSlength, Model) %>% 
  summarize(LE_pp=length(which(minci>0))/length(minci),
            LE_pp.01=length(which(minci>0.01))/length(minci),
            LEs_pp.05=length(which(minci>0.05))/length(minci))

#plots

#proportions, individual models
ggplot(sims_summary, aes(x=Model, y=LE_pp.01, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel2) + 
  geom_bar(stat = "identity", color="black") + classic + xlabvert + ylab("Proportion Classified Chaotic") +
  labs(fill="True\nModel\nDynamics") + 
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
ggsave("./figures/SimClassLE.pdf", width = 12, height = 7)

#different arrangement
ggplot(sims_summary, aes(x=NoiseLevel2, y=LE_pp.01, fill=Classification)) +
  facet_grid(TSlength~Classification) + 
  geom_bar(stat = "identity", position = position_dodge2(), show.legend = F) + theme_bw() + xlabvert


#plot values

#LE all models grouped together
ggplot(sims_d, aes(x=NoiseLevel2, y=minci, color=Classification)) +
  facet_grid(TSlength~., scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.05), alpha=0.4) + 
  theme_bw() + xlab("Noise Level") + legalpha
#violin plots
ggplot(sims_d, aes(x=factor(NoiseLevel2), y=minci, fill=Classification)) +
  facet_grid(TSlength~., scales = "free_y") + geom_hline(yintercept = 0) +
  geom_violin(position = position_dodge(0.9), scale = "width") + 
  theme_bw() + xlab("Noise Level")

#LE individual models
ggplot(sims_d, aes(x=Model, y=minci, color=Classification)) +
  facet_grid(TSlength~NoiseLevel2, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02), alpha=0.3) + theme_bw() + xlabvert

#R squared, E, theta for individual models
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
  group_by(Classification,LEsign,NoiseLevel2,TSlength) %>% 
  summarize(n=n()) %>% ungroup() %>% 
  complete(Classification, nesting(LEsign,NoiseLevel2,TSlength), 
           fill=list(n=0)) %>% 
  group_by(NoiseLevel2,TSlength, Classification) %>% 
  mutate(proportion=n/sum(n))
ggplot(sims_summary2, aes(x=LEsign, y=Classification, fill=proportion)) +
  facet_grid(TSlength~NoiseLevel2) + 
  geom_tile(stat = "identity") + classic +
  geom_text(aes(label=round(proportion,2)), color="white") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  ylab("True Dynamics") + xlab("Classification") +
  labs(fill="Proportion") 
#overall prop correct classification
sims_d %>% select(-data) %>% 
  mutate(Classification2=ifelse(Classification=="chaotic", "chaotic", "not chaotic")) %>% 
  group_by(Classification2,LEsign) %>% 
  summarize(n=n()) %>% ungroup() %>% 
  complete(Classification2, nesting(LEsign), 
           fill=list(n=0)) %>% 
  group_by(Classification2) %>% 
  mutate(proportion=n/sum(n))


#regression method
#class LEs
sims_d$LEregmin=sims_d$LEreg-1.96*sims_d$LEreg_se
sims_d$LEregsign=ifelse(sims_d$LEreg>0.01, "chaotic", "not chaotic")
sims_d$LEregsign2=ifelse(sims_d$LEregmin>0.01, "chaotic", "not chaotic")

#prop correct by model
sims_summaryreg=sims_d %>% select(-data) %>% 
  group_by(Classification,NoiseLevel2,TSlength, Model) %>% 
  summarize(LEreg_pp=length(which(LEregsign=="chaotic"))/length(LEregsign),
            LEreg_pp2=length(which(LEregsign2=="chaotic"))/length(LEregsign2))
ggplot(sims_summaryreg, aes(x=Model, y=LEreg_pp2, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel2) + 
  geom_bar(stat = "identity", color="black") + classic + xlabvert + ylab("Proportion Classified Chaotic") +
  labs(fill="True\nModel\nDynamics") + 
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))

#LE all models together
ggplot(sims_d, aes(x=NoiseLevel2, y=LEregmin, color=Classification)) +
  facet_grid(TSlength~.) + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.02)) + theme_bw()
#LE for each model
ggplot(sims_d, aes(x=Model, y=LEregmin, color=Classification)) +
  facet_grid(TSlength~NoiseLevel, scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point() + #geom_boxplot() + 
  theme_bw() + xlabvert

sims_summaryreg2=sims_d %>% select(-data) %>% 
  group_by(Classification,LEregsign2,NoiseLevel2,TSlength) %>% 
  summarize(n=n()) %>% ungroup() %>% 
  complete(Classification, nesting(LEregsign2,NoiseLevel2,TSlength), 
           fill=list(n=0)) %>% 
  group_by(NoiseLevel2,TSlength, Classification) %>% 
  mutate(proportion=n/sum(n))
ggplot(sims_summaryreg2, aes(x=LEregsign2, y=Classification, fill=proportion)) +
  facet_grid(TSlength~NoiseLevel2) + 
  geom_tile(stat = "identity") + classic +
  geom_text(aes(label=round(proportion,2)), color="white") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  ylab("True Dynamics") + xlab("Classification") +
  labs(fill="Proportion")

#overall prop correct classification
sims_d %>% select(-data) %>% 
  mutate(Classification2=ifelse(Classification=="chaotic", "chaotic", "not chaotic")) %>% 
  group_by(Classification2,LEregsign2) %>% 
  summarize(n=n()) %>% ungroup() %>% 
  complete(Classification2, nesting(LEregsign2), 
           fill=list(n=0)) %>% 
  group_by(Classification2) %>% 
  mutate(proportion=n/sum(n))

#validation data ####
#reclass noise level for stochastic ts
sims_v$NoiseLevel2=ifelse(sims_v$NoiseLevel==0, 0.01, sims_v$NoiseLevel)

#class LEs
sims_v$LEsign=ifelse(sims_v$minci>0.01, "chaotic", "not chaotic")

sims_summary=sims_v %>% select(-data) %>% 
  group_by(Classification,NoiseLevel2,TSlength, Model) %>% 
  summarize(LE_pp=length(which(minci>0))/length(minci),
            LE_pp.01=length(which(minci>0.01))/length(minci),
            LEs_pp.05=length(which(minci>0.05))/length(minci))

#plots

#proportions, individual models
ggplot(sims_summary, aes(x=Model, y=LE_pp.01, fill=Classification)) +
  facet_grid(TSlength~NoiseLevel2) + 
  geom_bar(stat = "identity", color="black") + classic + xlabvert + ylab("Proportion Classified Chaotic") +
  labs(fill="True\nModel\nDynamics") + 
  scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
