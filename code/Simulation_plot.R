# Plot simuation results

library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")
library(cowplot)

source("./code/ggplot themes rogers.R")

#load data ####
sims_d=read.csv("./data/sims_test_results_allmethods.csv")
sims_v=read.csv("./data/sims_validation_results_allmethods.csv")
sims_v2=read.csv("./data/sims_validation_results_forcedAR_allmethods.csv")
sims_v=rbind(sims_v, sims_v2)

#set model order
modelorder=unique(arrange(sims_d, Classification, Model)$Model)
sims_d$Model=factor(sims_d$Model, levels=modelorder)
modelorderv=unique(arrange(sims_v, Classification, Model)$Model)
sims_v$Model=factor(sims_v$Model, levels=modelorderv)

#convert to long format
sims_long=gather(sims_d, Method, Methodclass, LEregclass:DTclass)
sims_vlong=gather(sims_v, Method, Methodclass, LEregclass:DTclass)

#overall prop correct classification ####
sims_long %>% #test
  group_by(Method, Classification2, Methodclass) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method, nesting(Classification2, Methodclass), fill=list(n=0)) %>% 
  group_by(Method, Classification2) %>% mutate(proportion=n/sum(n)) %>% as.data.frame()
sims_vlong %>% #validation
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
  mutate(Method2=sub("class","", Method),Method2=recode(Method2, LE="JLE", LEreg="DLE"),Method2=factor(Method2, levels=c("DT","HVA","PE","RQA","JLE","DLE")))
mar2=sims_long %>%
  group_by(Method, NoiseLevel2) %>% summarize(errorrate=sum(error)/length(error)) %>% as.data.frame() %>% mutate(Dataset="Test Dataset") %>% 
  mutate(Method2=sub("class","", Method),Method2=recode(Method2, LE="JLE", LEreg="DLE"),Method2=factor(Method2, levels=c("DT","HVA","PE","RQA","JLE","DLE")))
mar1v=sims_vlong %>%
  group_by(Method, TSlength) %>% summarize(errorrate=sum(error)/length(error)) %>% as.data.frame() %>% mutate(Dataset="Validation Dataset") %>% 
  mutate(Method2=sub("class","", Method),Method2=recode(Method2, LE="JLE", LEreg="DLE"),Method2=factor(Method2, levels=c("DT","HVA","PE","RQA","JLE","DLE")))
mar2v=sims_vlong %>%
  group_by(Method, NoiseLevel2) %>% summarize(errorrate=sum(error)/length(error)) %>% as.data.frame() %>% mutate(Dataset="Validation Dataset") %>% 
  mutate(Method2=sub("class","", Method),Method2=recode(Method2, LE="JLE", LEreg="DLE"),Method2=factor(Method2, levels=c("DT","HVA","PE","RQA","JLE","DLE")))
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

#test data ####

#proportion correct classifications, across all models/methods ####
summary=sims_long %>% 
  group_by(Method,Classification,Methodclass,NoiseLevel2,TSlength) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method,nesting(Classification, Methodclass,NoiseLevel2,TSlength), fill=list(n=0)) %>% 
  group_by(NoiseLevel2,TSlength, Classification, Method) %>% 
  mutate(proportion=n/sum(n), Method2=sub("class","", Method), 
         Method2=recode(Method2, LE="JLE", LEreg="DLE"),
         Method2=factor(Method2, levels=c("DLE","JLE","RQA","PE","HVA","DT")))
ggplot(filter(summary, Methodclass=="chaotic"), aes(x=factor(NoiseLevel2), y=Classification, fill=proportion)) +
  facet_grid(Method2~TSlength) + geom_tile(stat = "identity") + 
  geom_text(aes(label=round(proportion,2)), color="white", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="True Dynamics", x="Noise Level", fill="Proportion Classified Chaotic", title = "Time Series Length") +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom",
                                panel.spacing.x = unit(0.2,"lines")) 
ggsave("./figures/SimResults.png", width = 7, height = 7)

#proportion correct classifications, individual models, top 3 methods ####
summary2=sims_long %>% 
  group_by(Method,Classification,NoiseLevel2,TSlength, Model) %>% 
  summarize(proportion=length(which(Methodclass=="chaotic"))/length(Methodclass)) %>% 
  mutate(Method2=sub("class","", Method), 
         Method2=recode(Method2, LE="JLE", LEreg="DLE"),
         Method2=factor(Method2, levels=c("DLE","JLE","RQA","PE","HVA","DT")))
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

#validation data ####

#proportion correct classifications, across all models/methods ####
summary=sims_vlong %>% 
  group_by(Method,Classification,Methodclass,NoiseLevel2,TSlength) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method,nesting(Classification, Methodclass,NoiseLevel2,TSlength), fill=list(n=0)) %>% 
  group_by(NoiseLevel2,TSlength, Classification, Method) %>% 
  mutate(proportion=n/sum(n), Method2=sub("class","", Method), 
         Method2=recode(Method2, LE="JLE", LEreg="DLE"),
         Method2=factor(Method2, levels=c("DLE", "JLE","RQA","PE","HVA","DT")))
ggplot(filter(summary, Methodclass=="chaotic"), aes(x=factor(NoiseLevel2), y=Classification, fill=proportion)) +
  facet_grid(Method2~TSlength) + geom_tile(stat = "identity") + 
  geom_text(aes(label=round(proportion,2)), color="white", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="True Dynamics", x="Noise Level", fill="Proportion Classified Chaotic", title = "Time Series Length") +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom",
                                panel.spacing.x = unit(0.2,"lines")) 
ggsave("./figures/SimResults_validation.png", width = 7, height = 7)

#proportion correct classifications, individual models, top 3 methods ####
summary2=sims_vlong %>% 
  group_by(Method,Classification,NoiseLevel2,TSlength, Model) %>% 
  summarize(proportion=length(which(Methodclass=="chaotic"))/length(Methodclass)) %>% 
  mutate(Method2=sub("class","", Method), 
         Method2=recode(Method2, LE="JLE", LEreg="DLE"),
         Method2=factor(Method2, levels=c("DLE", "JLE","RQA","PE","HVA","DT")))
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



# misc other plots (not used) ####

plotprop1=function(method, title) {
  ggplot(filter(summary, Method==method), aes(x=Methodclass, y=Classification, fill=proportion)) +
    facet_grid(TSlength~NoiseLevel2) + geom_tile(stat = "identity") + 
    geom_text(aes(label=round(proportion,2)), color="white") +
    classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
    labs(y="True Dynamics", x="Classification", fill="Proportion", title = title) 
}
plotprop1("LEregclass", "Regression LE")
plotprop1("LEclass", "JLE")
plotprop1("RQAclass", "RQA")
plotprop1("PEclass", "PE")
plotprop1("HVAclass", "HVA")
plotprop1("DTclass", "DT")

plotprop2=function(method, title) {
  ggplot(filter(summary2, Method==method), aes(x=Model, y=LE_pp, fill=Classification)) +
    facet_grid(TSlength~NoiseLevel2) + classic + xlabvert +
    geom_bar(stat = "identity", color="black") + 
    labs(fill="True\nModel\nDynamics", y="Proportion Classified Chaotic", title=title) + 
    scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
}
plotprop2("LEregclass", "Regression LE")
plotprop2("LEclass", "JLE")
plotprop2("RQAclass", "RQA")
plotprop2("PEclass", "PE")
plotprop2("HVAclass", "HVA")
plotprop2("DTclass", "DT")

#test and validation together
#proportion correct classifications, across all models
summary=rbind(sims_long, sims_vlong) %>% 
  group_by(Method,Classification,Methodclass,NoiseLevel2,TSlength) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method,nesting(Classification, Methodclass,NoiseLevel2,TSlength), fill=list(n=0)) %>% 
  group_by(NoiseLevel2,TSlength, Classification, Method) %>% mutate(proportion=n/sum(n))
plotprop1("LEregclass", "Regression LE")
plotprop1("LEclass", "JLE")
plotprop1("RQAclass", "RQA")
plotprop1("PEclass", "PE")
plotprop1("HVAclass", "HVA")
plotprop1("DTclass", "DT")

#proportion correct classifications, individual models
summary2=rbind(sims_long, sims_vlong) %>% 
  group_by(Method,Classification,NoiseLevel2,TSlength, Model) %>% 
  summarize(LE_pp=length(which(Methodclass=="chaotic"))/length(Methodclass))
modelorder2=as.character(unique(arrange(summary2, Classification, Model)$Model))
summary2$Model=factor(summary2$Model, levels=modelorder2)
plotprop2("LEregclass", "Regression LE")
plotprop2("LEclass", "JLE")
plotprop2("RQAclass", "RQA")
plotprop2("PEclass", "PE")
plotprop2("HVAclass", "HVA")
plotprop2("DTclass", "DT")

#different arrangement
ggplot(sims_summary2, aes(x=NoiseLevel2, y=LE_pp, fill=Classification)) +
  facet_grid(TSlength~Classification) + 
  geom_bar(stat = "identity", position = position_dodge2(), show.legend = F) + theme_bw() + xlabvert

#LE across all models
ggplot(sims_d, aes(x=NoiseLevel2, y=LEmin, color=Classification)) +
  facet_grid(TSlength~., scales = "free_y") + geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.05), alpha=0.4) + 
  theme_bw() + xlab("Noise Level") + legalpha
#violin plots
ggplot(sims_d, aes(x=factor(NoiseLevel2), y=LEmin, fill=Classification)) +
  facet_grid(TSlength~., scales = "free_y") + geom_hline(yintercept = 0) +
  geom_violin(position = position_dodge(0.9), scale = "width") + 
  theme_bw() + xlab("Noise Level")

#LE for individual models
ggplot(sims_d, aes(x=Model, y=LEmin, color=Classification)) +
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
  print(ggplot(filter(sims_d, TSlength==tslengths[i]), aes(x=LEmin, fill=Classification)) +
          facet_grid(Classification~NoiseLevel2, scales = "free_y") + geom_hline(yintercept = 0) +
          geom_histogram(boundary = 0.01, binwidth = 0.1, show.legend = F) + geom_vline(xintercept = 0) +
          theme_bw() + ggtitle(paste("TSlength =", tslengths[i])) + xlim(c(-1,1)))
}

