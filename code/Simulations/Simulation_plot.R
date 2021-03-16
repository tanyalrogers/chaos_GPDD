# Plots of simuation results
# Tanya Rogers

library(ggplot2)
library(dplyr)
library(tidyr)
library(rEDM)
library(cowplot)

source("./code/Simulations/ggplot_themes_rogers.R")

#load data
sims_d=read.csv("./data/sims_test_results.csv")
sims_v=read.csv("./data/sims_validation_results.csv")
sims_n=read.csv("./data/sims_noise_results.csv")
sims_d2=read.csv("./data/sims_test_results_othermethods.csv")
sims_v2=read.csv("./data/sims_validation_results_othermethods.csv")
sims_n2=read.csv("./data/sims_noise_results_othermethods.csv")

sims_d=left_join(sims_d, sims_d2)
sims_v=left_join(sims_v, sims_v2)
sims_n=left_join(sims_n, sims_n2)

#specify dynamics of noise models
sims_n=unite(sims_n,Model,Model,Classification, remove = F)

#set model order
modelorder=unique(arrange(sims_d, Classification, Model)$Model)
sims_d$Model=factor(sims_d$Model, levels=modelorder)
modelorderv=unique(arrange(sims_v, Classification, Model)$Model)
sims_v$Model=factor(sims_v$Model, levels=modelorderv)
modelordern=unique(arrange(sims_n, Classification, Model)$Model)
sims_n$Model=factor(sims_n$Model, levels=modelordern)

#convert to long format
sims_long=gather(sims_d, Method, Methodclass, LEregclass:DTclass)
sims_vlong=gather(sims_v, Method, Methodclass, LEregclass:DTclass)
sims_nlong=gather(sims_n, Method, Methodclass, LEregclass:DTclass)

#overall prop correct classification ####
sims_long %>% #test
  group_by(Method, Classification2, Methodclass) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method, nesting(Classification2, Methodclass), fill=list(n=0)) %>% 
  group_by(Method, Classification2) %>% mutate(proportion=n/sum(n)) %>% as.data.frame()
sims_vlong %>% #validation
  group_by(Method, Classification2, Methodclass) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method, nesting(Classification2, Methodclass), fill=list(n=0)) %>% 
  group_by(Method, Classification2) %>% mutate(proportion=n/sum(n)) %>% as.data.frame()
sims_nlong %>% #noise
  group_by(Method, Classification2, Methodclass) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method, nesting(Classification2, Methodclass), fill=list(n=0)) %>% 
  group_by(Method, Classification2) %>% mutate(proportion=n/sum(n)) %>% as.data.frame()
#test and validation together
sims_long %>% rbind(sims_vlong) %>% 
  group_by(Method, Classification2, Methodclass) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method, nesting(Classification2, Methodclass), fill=list(n=0)) %>% 
  group_by(Method, Classification2) %>% mutate(proportion=n/sum(n)) %>% as.data.frame()

#error rate, marginalized by tslength, noiselevel, test and validation ####
sims_long$error=ifelse(sims_long$Classification2==sims_long$Methodclass,0,1)
sims_vlong$error=ifelse(sims_vlong$Classification2==sims_vlong$Methodclass,0,1)

mar1=sims_long %>%
  group_by(Method, TSlength) %>% summarize(errorrate=sum(error)/length(error)) %>% as.data.frame() %>% mutate(Dataset="Test Dataset") %>% 
  mutate(Method2=sub("class","", Method),Method2=recode(Method2, LE="JLE", LEreg="DLE", DT="CDT", HVA="HVG"),Method2=factor(Method2, levels=c("CDT","HVG","PE","RQA","JLE","DLE")))
mar2=sims_long %>%
  group_by(Method, NoiseLevel2) %>% summarize(errorrate=sum(error)/length(error)) %>% as.data.frame() %>% mutate(Dataset="Test Dataset") %>% 
  mutate(Method2=sub("class","", Method),Method2=recode(Method2, LE="JLE", LEreg="DLE", DT="CDT", HVA="HVG"),Method2=factor(Method2, levels=c("CDT","HVG","PE","RQA","JLE","DLE")))
mar1v=sims_vlong %>%
  group_by(Method, TSlength) %>% summarize(errorrate=sum(error)/length(error)) %>% as.data.frame() %>% mutate(Dataset="Validation Dataset") %>% 
  mutate(Method2=sub("class","", Method),Method2=recode(Method2, LE="JLE", LEreg="DLE", DT="CDT", HVA="HVG"),Method2=factor(Method2, levels=c("CDT","HVG","PE","RQA","JLE","DLE")))
mar2v=sims_vlong %>%
  group_by(Method, NoiseLevel2) %>% summarize(errorrate=sum(error)/length(error)) %>% as.data.frame() %>% mutate(Dataset="Validation Dataset") %>% 
  mutate(Method2=sub("class","", Method),Method2=recode(Method2, LE="JLE", LEreg="DLE", DT="CDT", HVA="HVG"),Method2=factor(Method2, levels=c("CDT","HVG","PE","RQA","JLE","DLE")))
mar1=rbind(mar1, mar1v)
mar2=rbind(mar2, mar2v)

er1=ggplot(mar1, aes(x=factor(TSlength), y=Method2, fill=errorrate)) +
  facet_grid(.~Dataset) + 
  geom_tile(stat = "identity") + 
  geom_text(aes(label=round(errorrate,2)), color="black", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="Method", x="Time Series Length", fill="Classification Error Rate") +
  scale_fill_distiller(palette = "Reds", direction = 1, limits=range(mar1$errorrate,mar2$errorrate)) + 
  removefacetbackground + 
  theme(plot.title = element_text(hjust = 0.5, size=11), 
        strip.text = element_text(size=10), 
        legend.position = "bottom", panel.spacing.x = unit(0.75,"lines")) 
er2=ggplot(mar2, aes(x=factor(NoiseLevel2), y=Method2, fill=errorrate)) +
  facet_grid(.~Dataset) + 
  geom_tile(stat = "identity") + 
  geom_text(aes(label=round(errorrate,2)), color="black", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="Method", x="Noise Level", fill="Classification Error Rate") +
  scale_fill_distiller(palette = "Reds", direction = 1, limits=range(mar1$errorrate,mar2$errorrate)) + 
  removefacetbackground + 
  theme(plot.title = element_text(hjust = 0.5, size=11), 
        strip.text = element_text(size=10), 
        legend.position = "bottom", panel.spacing.x = unit(0.75,"lines")) 

legend <- get_legend(er1)
plot_grid(er1+ theme(legend.position="none"),er2+ theme(legend.position="none"), legend, nrow = 3, rel_heights = c(1,1,0.2))
ggsave("./figures/marginal_errorrates.png", width = 6, height = 6)

#by ts length and error rate (not used)
mar3=sims_long %>%
  group_by(Method, TSlength, NoiseLevel2) %>% summarize(errorrate=sum(error)/length(error)) %>% as.data.frame() %>% 
  mutate(Method2=sub("class","", Method), Method2=recode(Method2, LE="JLE", LEreg="DLE"),Method2=factor(Method2, levels=c("DLE", "JLE","RQA","PE","HVA","DT")))
ggplot(mar3, aes(x=factor(NoiseLevel2), y=factor(TSlength), fill=errorrate)) +
  facet_wrap(.~Method2, ncol=3) + 
  geom_tile(stat = "identity") + 
  geom_text(aes(label=round(errorrate,2)), color="black", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="Time Series Length", x="Noise Level", fill="Classification Error Rate") +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom", panel.spacing = unit(0.5,"lines")) 
mar3v=sims_vlong %>%
  group_by(Method, TSlength, NoiseLevel2) %>% summarize(errorrate=sum(error)/length(error)) %>% as.data.frame() %>% 
  mutate(Method2=sub("class","", Method), Method2=recode(Method2, LE="JLE", LEreg="DLE"),Method2=factor(Method2, levels=c("DLE", "JLE","RQA","PE","HVA","DT")))
ggplot(mar3v, aes(x=factor(NoiseLevel2), y=factor(TSlength), fill=errorrate)) +
  facet_wrap(.~Method2, ncol=3) + 
  geom_tile(stat = "identity") + 
  geom_text(aes(label=round(errorrate,2)), color="black", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="Time Series Length", x="Noise Level", fill="Classification Error Rate") +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom", panel.spacing = unit(0.5,"lines")) 

#Test dataset ####

#proportion correct classifications, across all models/methods
summary=sims_long %>% 
  group_by(Method,Classification,Methodclass,NoiseLevel2,TSlength) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method,nesting(Classification, Methodclass,NoiseLevel2,TSlength), fill=list(n=0)) %>% 
  group_by(NoiseLevel2,TSlength, Classification, Method) %>% 
  mutate(proportion=n/sum(n), Method2=sub("class","", Method), 
         Method2=recode(Method2, LE="JLE", LEreg="DLE", DT="CDT", HVA="HVG"),
         Method2=factor(Method2, levels=c("DLE","JLE","RQA","PE","HVG","CDT")))
ggplot(filter(summary, Methodclass=="chaotic"), aes(x=factor(NoiseLevel2), y=Classification, fill=proportion)) +
  facet_grid(Method2~TSlength) + geom_tile(stat = "identity") + 
  geom_text(aes(label=round(proportion,2)), color="white", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="True Dynamics", x="Noise Level", fill="Proportion Classified Chaotic", title = "Time Series Length") +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom",
                                panel.spacing.x = unit(0.2,"lines")) 
ggsave("./figures/SimResults.png", width = 7, height = 7)

#proportion correct classifications, individual models, top 3 methods
summary2=sims_long %>% 
  group_by(Method,Classification,NoiseLevel2,TSlength, Model) %>% 
  summarize(proportion=length(which(Methodclass=="chaotic"))/length(Methodclass)) %>% 
  mutate(Method2=sub("class","", Method), 
         Method2=recode(Method2, LE="JLE", LEreg="DLE", DT="CDT", HVA="HVG"),
         Method2=factor(Method2, levels=c("DLE","JLE","RQA","PE","HVG","CDT")))
ggplot(filter(summary2, Method2 %in% c("JLE", "RQA", "PE")), aes(x=factor(NoiseLevel2), y=Model, fill=proportion)) +
  facet_grid(Method2~TSlength) + geom_tile(stat = "identity") + 
  geom_hline(yintercept = c(7.5,14.5), color="white") +
  #geom_text(aes(label=round(proportion,2)), color="white", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="Model", x="Noise Level", fill="Proportion Classified Chaotic", title = "Time Series Length") +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom",
                                axis.text.y = element_text(size=8, color=rep(c("black","royalblue","seagreen"), each=7)),
                                panel.spacing.x = unit(0.2,"lines")) 
ggsave("./figures/SimResultsIndivModels.png", width = 7, height = 8)

#Validation data ####

#proportion correct classifications, across all models/methods
summary=sims_vlong %>% 
  group_by(Method,Classification,Methodclass,NoiseLevel2,TSlength) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method,nesting(Classification, Methodclass,NoiseLevel2,TSlength), fill=list(n=0)) %>% 
  group_by(NoiseLevel2,TSlength, Classification, Method) %>% 
  mutate(proportion=n/sum(n), Method2=sub("class","", Method), 
         Method2=recode(Method2, LE="JLE", LEreg="DLE", DT="CDT", HVA="HVG"),
         Method2=factor(Method2, levels=c("DLE", "JLE","RQA","PE","HVG","CDT")))
ggplot(filter(summary, Methodclass=="chaotic"), aes(x=factor(NoiseLevel2), y=Classification, fill=proportion)) +
  facet_grid(Method2~TSlength) + geom_tile(stat = "identity") + 
  geom_text(aes(label=round(proportion,2)), color="white", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="True Dynamics", x="Noise Level", fill="Proportion Classified Chaotic", title = "Time Series Length") +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom",
                                panel.spacing.x = unit(0.2,"lines")) 
ggsave("./figures/SimResults_validation.png", width = 7, height = 7)

#proportion correct classifications, individual models, top 3 methods
summary2=sims_vlong %>% 
  group_by(Method,Classification,NoiseLevel2,TSlength, Model) %>% 
  summarize(proportion=length(which(Methodclass=="chaotic"))/length(Methodclass)) %>% 
  mutate(Method2=sub("class","", Method), 
         Method2=recode(Method2, LE="JLE", LEreg="DLE", DT="CDT", HVA="HVG"),
         Method2=factor(Method2, levels=c("DLE", "JLE","RQA","PE","HVG","CDT")))
ggplot(filter(summary2, Method2 %in% c("JLE", "RQA", "PE")), aes(x=factor(NoiseLevel2), y=Model, fill=proportion)) +
  facet_grid(Method2~TSlength) + geom_tile(stat = "identity") + 
  geom_hline(yintercept = c(3.5,6.5), color="white") +
  #geom_text(aes(label=round(proportion,2)), color="white", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="Model", x="Noise Level", fill="Proportion Classified Chaotic", title = "Time Series Length") +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom",
                                axis.text.y = element_text(size=8, color=c(rep(c("black","royalblue","seagreen"), each=3),"seagreen")),
                                panel.spacing.x = unit(0.2,"lines")) 
ggsave("./figures/SimResultsIndivModels_validation.png", width = 7, height = 6)

#Observation and process noise dataset ####

#proportion correct classifications, individual models, top 3 methods
summary2=sims_nlong %>% 
  group_by(Method,Classification,ObsNoise,ProcessNoise, Model) %>% 
  summarize(proportion=length(which(Methodclass=="chaotic"))/length(Methodclass)) %>% 
  mutate(Method2=sub("class","", Method),
         Method2=recode(Method2, LE="JLE", LEreg="DLE", DT="CDT", HVA="HVG"),
         Method2=factor(Method2, levels=c("DLE", "JLE","RQA","PE","HVG","CDT")))
ggplot(filter(summary2, Method2 %in% c("JLE", "RQA", "PE")), aes(x=factor(ObsNoise), y=Model, fill=proportion)) +
  facet_grid(Method2~ProcessNoise) + geom_tile(stat = "identity") + 
  geom_hline(yintercept = 3.5, color="white") +
  #geom_text(aes(label=round(proportion,2)), color="white", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="Model", x="Observation Noise Level", fill="Proportion Classified Chaotic", title = "Process Noise Level") +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom",
                                axis.text.y = element_text(size=8, color=c(rep(c("black","royalblue"), each=3))),
                                panel.spacing.x = unit(0.2,"lines")) 
ggsave("./figures/SimResultsIndivModels_noise.png", width = 7, height = 5)
