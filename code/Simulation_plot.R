# Plot simuation results

library("ggplot2")
library("dplyr")
library("tidyr")
library("rEDM")

source("./code/ggplot themes rogers.R")

#test data ####

#load results
#load("./data/sims_results_update.Rdata")
sims_d=read.csv("./data/sims_test_results.csv", stringsAsFactors = F)
modelorder=unique(arrange(sims_d, Classification, Model)$Model)
sims_d$Model=factor(sims_d$Model, levels=modelorder)
#reclass noise level for stochastic ts
sims_d$NoiseLevel2=ifelse(sims_d$NoiseLevel==0, 0.01, sims_d$NoiseLevel)
sims_d$Classification2=ifelse(sims_d$Classification=="chaotic", "chaotic", "not chaotic")

#class LEs regression method
sims_d$LEregmin=sims_d$LEreg-1.96*sims_d$LEreg_se
sims_d$LEregclass=ifelse(sims_d$LEregmin>0.01, "chaotic", "not chaotic")
#class LEs Jacobian method
sims_d$LEclass=ifelse(sims_d$LEmin>0.01, "chaotic", "not chaotic")
#class other methods
RQAclass=read.csv("./data/RQAclassification.csv") %>% gather(SimNumber, RQAclass, Sim.1:Sim.100)
PEclass=read.csv("./data/PEclassification.csv") %>% gather(SimNumber, PEclass, Sim.1:Sim.100)
HVAclass=read.csv("./data/HVAclassification.csv") %>% gather(SimNumber, HVAclass, Sim.1:Sim.100)
DTclass=read.csv("./data/DTclassification.csv") %>% gather(SimNumber, DTclass, Sim.1:Sim.100)
#join to main table
sims_d=left_join(sims_d, RQAclass) %>% left_join(PEclass) %>% 
  left_join(HVAclass) %>% left_join(DTclass)
#long format
sims_long=gather(sims_d, Method, Methodclass, LEregclass:DTclass)

#overall prop correct classification
sims_long %>%
  group_by(Method, Classification2, Methodclass) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method, nesting(Classification2, Methodclass), fill=list(n=0)) %>% 
  group_by(Method, Classification2) %>% mutate(proportion=n/sum(n)) %>% as.data.frame()

#plots

#error rate, marginalized by tslength, noiselevel
sims_long$error=ifelse(sims_long$Classification2==sims_long$Methodclass,0,1)
mar1=sims_long %>%
  group_by(Method, TSlength) %>% summarize(errorrate=sum(error)/length(error)) %>% as.data.frame() %>% 
  mutate(Method2=sub("class","", Method),Method2=recode(Method2, LE="Jacobian LE", LEreg="Direct LE"),Method2=factor(Method2, levels=c("DT","HVA","PE","RQA","Jacobian LE","Direct LE")))
mar2=sims_long %>%
  group_by(Method, NoiseLevel2) %>% summarize(errorrate=sum(error)/length(error)) %>% as.data.frame() %>% 
  mutate(Method2=sub("class","", Method),Method2=recode(Method2, LE="Jacobian LE", LEreg="Direct LE"),Method2=factor(Method2, levels=c("DT","HVA","PE","RQA","Jacobian LE","Direct LE")))

#barplot
ggplot(mar1, aes(x=factor(TSlength), y=errorrate)) +
  facet_wrap(.~Method2, ncol=3) + 
  geom_bar(stat = "identity", fill="firebrick") + 
  classic + scale_y_continuous(expand = c(0,0), limits = c(0,0.55)) +
  labs(x="Time Series Length", y="Classification Error Rate") + removefacetbackground
 
ggplot(mar1, aes(x=factor(TSlength), y=Method2, fill=errorrate)) +
  #facet_grid(Method2~TSlength) + 
  geom_tile(stat = "identity") + 
  geom_text(aes(label=round(errorrate,2)), color="black", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="Method", x="Time Series Length", fill="Classification Error Rate", title = "Test Dataset") +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom", panel.spacing.x = unit(0.2,"lines")) 
ggplot(mar2, aes(x=factor(NoiseLevel2), y=Method2, fill=errorrate)) +
  #facet_grid(Method2~TSlength) + 
  geom_tile(stat = "identity") + 
  geom_text(aes(label=round(errorrate,2)), color="black", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="Method", x="Noise Level", fill="Classification Error Rate", title = "Test Dataset") +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom", panel.spacing.x = unit(0.2,"lines")) 

mar3=sims_long %>%
  group_by(Method, TSlength, NoiseLevel2) %>% summarize(errorrate=sum(error)/length(error)) %>% as.data.frame() %>% 
  mutate(Method2=sub("class","", Method), Method2=recode(Method2, LE="Jacobian LE", LEreg="Direct LE"),Method2=factor(Method2, levels=c("Direct LE", "Jacobian LE","RQA","PE","HVA","DT")))
ggplot(mar3, aes(x=factor(NoiseLevel2), y=factor(TSlength), fill=errorrate)) +
  facet_wrap(.~Method2, ncol=3) + 
  geom_tile(stat = "identity") + 
  geom_text(aes(label=round(errorrate,2)), color="black", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="Time Series Length", x="Noise Level", fill="Classification Error Rate") +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom", panel.spacing = unit(0.5,"lines")) 


#proportion correct classifications, across all models
summary=sims_long %>% 
  group_by(Method,Classification,Methodclass,NoiseLevel2,TSlength) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method,nesting(Classification, Methodclass,NoiseLevel2,TSlength), fill=list(n=0)) %>% 
  group_by(NoiseLevel2,TSlength, Classification, Method) %>% 
  mutate(proportion=n/sum(n), Method2=sub("class","", Method), 
         Method2=recode(Method2, LE="Jacobian LE", LEreg="Direct LE"),
         Method2=factor(Method2, levels=c("Direct LE", "Jacobian LE","RQA","PE","HVA","DT")))
plotprop1=function(method, title) {
  ggplot(filter(summary, Method==method), aes(x=Methodclass, y=Classification, fill=proportion)) +
  facet_grid(TSlength~NoiseLevel2) + geom_tile(stat = "identity") + 
  geom_text(aes(label=round(proportion,2)), color="white") +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="True Dynamics", x="Classification", fill="Proportion", title = title) 
}
plotprop1("LEregclass", "Regression LE")
plotprop1("LEclass", "Jacobian LE")
plotprop1("RQAclass", "RQA")
plotprop1("PEclass", "PE")
plotprop1("HVAclass", "HVA")
plotprop1("DTclass", "DT")

#all methods together ####
ggplot(filter(summary, Methodclass=="chaotic"), aes(x=factor(NoiseLevel2), y=Classification, fill=proportion)) +
  facet_grid(Method2~TSlength) + geom_tile(stat = "identity") + 
  geom_text(aes(label=round(proportion,2)), color="white", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="True Dynamics", x="Noise Level", fill="Proportion Classified Chaotic", title = "Time Series Length") +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom",
                                panel.spacing.x = unit(0.2,"lines")) 
ggsave("./figures/SimResults.png", width = 7, height = 7)

#proportion correct classifications, individual models
summary2=sims_long %>% 
  group_by(Method,Classification,NoiseLevel2,TSlength, Model) %>% 
  summarize(proportion=length(which(Methodclass=="chaotic"))/length(Methodclass)) %>% 
  mutate(Method2=sub("class","", Method), 
         Method2=recode(Method2, LE="Jacobian LE", LEreg="Direct LE"),
         Method2=factor(Method2, levels=c("Direct LE", "Jacobian LE","RQA","PE","HVA","DT")))
plotprop2=function(method, title) {
  ggplot(filter(summary2, Method==method), aes(x=Model, y=LE_pp, fill=Classification)) +
    facet_grid(TSlength~NoiseLevel2) + classic + xlabvert +
    geom_bar(stat = "identity", color="black") + 
    labs(fill="True\nModel\nDynamics", y="Proportion Classified Chaotic", title=title) + 
    scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0))
}
plotprop2("LEregclass", "Regression LE")
plotprop2("LEclass", "Jacobian LE")
plotprop2("RQAclass", "RQA")
plotprop2("PEclass", "PE")
plotprop2("HVAclass", "HVA")
plotprop2("DTclass", "DT")

#individual models ####
ggplot(filter(summary2, Method2 %in% c("Jacobian LE", "RQA", "PE")), aes(x=factor(NoiseLevel2), y=Model, fill=proportion)) +
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

sims_v=read.csv("./data/sims_validation_results.csv", stringsAsFactors = F)
modelorderv=unique(arrange(sims_v, Classification, Model)$Model)
sims_v$Model=factor(sims_v$Model, levels=modelorderv)
#reclass noise level for stochastic ts
sims_v$NoiseLevel2=ifelse(sims_v$NoiseLevel==0, 0.01, sims_v$NoiseLevel)
sims_v$Classification2=ifelse(sims_v$Classification=="chaotic", "chaotic", "not chaotic")

#class LEs regression method
sims_v$LEregmin=sims_v$LEreg-1.96*sims_v$LEreg_se
sims_v$LEregclass=ifelse(sims_v$LEregmin>0.01, "chaotic", "not chaotic")
#class LEs jacobian method
sims_v$LEclass=ifelse(sims_v$LEmin>0.01, "chaotic", "not chaotic")
#class other methods
RQAclassv=read.csv("./data/RQAclassification_validation.csv") %>% gather(SimNumber, RQAclass, Sim.1:Sim.100)
PEclassv=read.csv("./data/PEclassification_validation.csv") %>% gather(SimNumber, PEclass, Sim.1:Sim.100)
HVAclassv=read.csv("./data/HVAclassification_validation.csv") %>% gather(SimNumber, HVAclass, Sim.1:Sim.100)
DTclassv=read.csv("./data/DTclassification_validation.csv") %>% gather(SimNumber, DTclass, Sim.1:Sim.100)
#join to main table
sims_v=left_join(sims_v, RQAclassv) %>% left_join(PEclassv) %>%
  left_join(HVAclassv) %>% left_join(DTclassv)
#long format
sims_vlong=gather(sims_v, Method, Methodclass, LEregclass:DTclass)

#overall prop correct classification
sims_vlong %>%
  group_by(Method, Classification2, Methodclass) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method, nesting(Classification2, Methodclass), fill=list(n=0)) %>% 
  group_by(Method, Classification2) %>% mutate(proportion=n/sum(n)) %>% as.data.frame()

#plots
sims_vlong$error=ifelse(sims_vlong$Classification2==sims_vlong$Methodclass,0,1)
mar3v=sims_vlong %>%
  group_by(Method, TSlength, NoiseLevel2) %>% summarize(errorrate=sum(error)/length(error)) %>% as.data.frame() %>% 
  mutate(Method2=sub("class","", Method), Method2=recode(Method2, LE="Jacobian LE", LEreg="Direct LE"),Method2=factor(Method2, levels=c("Direct LE", "Jacobian LE","RQA","PE","HVA","DT")))
ggplot(mar3v, aes(x=factor(NoiseLevel2), y=factor(TSlength), fill=errorrate)) +
  facet_wrap(.~Method2, ncol=3) + 
  geom_tile(stat = "identity") + 
  geom_text(aes(label=round(errorrate,2)), color="black", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="Time Series Length", x="Noise Level", fill="Classification Error Rate") +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom", panel.spacing = unit(0.5,"lines")) 

#proportion correct classifications, across all models
summary=sims_vlong %>% 
  group_by(Method,Classification,Methodclass,NoiseLevel2,TSlength) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method,nesting(Classification, Methodclass,NoiseLevel2,TSlength), fill=list(n=0)) %>% 
  group_by(NoiseLevel2,TSlength, Classification, Method) %>% 
  mutate(proportion=n/sum(n), Method2=sub("class","", Method), 
         Method2=recode(Method2, LE="Jacobian LE", LEreg="Direct LE"),
         Method2=factor(Method2, levels=c("Direct LE", "Jacobian LE","RQA","PE","HVA","DT")))
plotprop1("LEregclass", "Regression LE")
plotprop1("LEclass", "Jacobian LE")
plotprop1("RQAclass", "RQA")
plotprop1("PEclass", "PE")
plotprop1("HVAclass", "HVA")
plotprop1("DTclass", "DT")

#all methods together validation ####
ggplot(filter(summary, Methodclass=="chaotic"), aes(x=factor(NoiseLevel2), y=Classification, fill=proportion)) +
  facet_grid(Method2~TSlength) + geom_tile(stat = "identity") + 
  geom_text(aes(label=round(proportion,2)), color="white", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="True Dynamics", x="Noise Level", fill="Proportion Classified Chaotic", title = "Time Series Length") +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom",
                                panel.spacing.x = unit(0.2,"lines")) 
ggsave("./figures/SimResults_validation.png", width = 7, height = 7)


#proportion correct classifications, individual models
summary2=sims_vlong %>% 
  group_by(Method,Classification,NoiseLevel2,TSlength, Model) %>% 
  summarize(proportion=length(which(Methodclass=="chaotic"))/length(Methodclass)) %>% 
  mutate(Method2=sub("class","", Method), 
         Method2=recode(Method2, LE="Jacobian LE", LEreg="Direct LE"),
         Method2=factor(Method2, levels=c("Direct LE", "Jacobian LE","RQA","PE","HVA","DT")))
plotprop2("LEregclass", "Regression LE")
plotprop2("LEclass", "Jacobian LE")
plotprop2("RQAclass", "RQA")
plotprop2("PEclass", "PE")
plotprop2("HVAclass", "HVA")
plotprop2("DTclass", "DT")

#individual models validation ####
ggplot(filter(summary2, Method2 %in% c("Jacobian LE", "RQA", "PE")), aes(x=factor(NoiseLevel2), y=Model, fill=proportion)) +
  facet_grid(Method2~TSlength) + geom_tile(stat = "identity") + 
  geom_hline(yintercept = c(3.5,6.5), color="white") +
  #geom_text(aes(label=round(proportion,2)), color="white", size=3) +
  classic + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="Model", x="Noise Level", fill="Proportion Classified Chaotic", title = "Time Series Length") +
  removefacetbackground + theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom",
                                axis.text.y = element_text(size=8, color=rep(c("black","royalblue","seagreen"), each=3)),
                                panel.spacing.x = unit(0.2,"lines")) 
ggsave("./figures/SimResultsIndivModels_validation.png", width = 7, height = 6)


#test and validation together ####

#overall prop correct classification
sims_long %>% rbind(sims_vlong) %>% 
  group_by(Method, Classification2, Methodclass) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method, nesting(Classification2, Methodclass), fill=list(n=0)) %>% 
  group_by(Method, Classification2) %>% mutate(proportion=n/sum(n)) %>% as.data.frame()

#proportion correct classifications, across all models
summary=rbind(sims_long, sims_vlong) %>% 
  group_by(Method,Classification,Methodclass,NoiseLevel2,TSlength) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method,nesting(Classification, Methodclass,NoiseLevel2,TSlength), fill=list(n=0)) %>% 
  group_by(NoiseLevel2,TSlength, Classification, Method) %>% mutate(proportion=n/sum(n))
plotprop1("LEregclass", "Regression LE")
plotprop1("LEclass", "Jacobian LE")
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
plotprop2("LEclass", "Jacobian LE")
plotprop2("RQAclass", "RQA")
plotprop2("PEclass", "PE")
plotprop2("HVAclass", "HVA")
plotprop2("DTclass", "DT")

# misc other plots ####

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

