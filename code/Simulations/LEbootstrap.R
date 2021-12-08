#LE bootstrapping
library(dplyr)
library(tidyr)
library(rEDM)
library(purrr)
library(furrr)
library(ggplot2)

source("./code/Methods/LEbootstrapfunctions.R")
source("./code/Methods/LE_ChaosDetectionMethods.R")

load(file = "./data/sims_results_update2.Rdata")

#subset for testing

SIMS=paste0("Sim.",1:20)

sims_sub=filter(sims_d, SimNumber %in% SIMS)
sims_results_sub=filter(sims_results, ID %in% unique(sims_sub$ID) & SimNumber %in% SIMS)
modelresultsbest=sims_results_sub$modelresultsbest
#test3=map(modelresultsbest,LEboot)

#original results (all) ####
#FPR
length(which(sims_sub$minci[sims_sub$Classification!="chaotic"]>0.01))/length(sims_sub$minci[sims_sub$Classification!="chaotic"])
#FNR
length(which(sims_sub$minci[sims_sub$Classification=="chaotic"]<0.01))/length(sims_sub$minci[sims_sub$Classification=="chaotic"])

tab2=data.frame(TSlength=sims_sub$TSlength,NoiseLevel=sims_sub$NoiseLevel, Classification=sims_sub$Classification, LEtest=sims_d$minci)
err=aggregate(LEtest~NoiseLevel*TSlength*Classification, data=tab2, FUN=function(x) {length(which(x>0.01))/length(x)})
err$NoiseLevel2=ifelse(err$NoiseLevel==0, 0.01, err$NoiseLevel)

ggplot(err, aes(x=factor(NoiseLevel2), y=factor(TSlength), fill=LEtest)) +
  facet_grid(.~Classification) +
  geom_tile() + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  geom_text(aes(label=round(LEtest,2)), color="black", size=3) +
  labs(y="Time Series Length", x="Noise Level", fill="Prop chaotic") +
  scale_fill_distiller(palette = "Blues", direction = 1, limits=c(0,1))

# GLE point estimate ####
LEpoint=map2_dbl(modelresultsbest,sims_results_sub$jacobians,LEglobal)
# save(LEpoint, file = "data/LEpoint.RData")

#FPR
length(which(LEpoint[sims_sub$Classification!="chaotic"]>0))/length(LEpoint[sims_sub$Classification!="chaotic"])
#FNR
length(which(LEpoint[sims_sub$Classification=="chaotic"]<0))/length(LEpoint[sims_sub$Classification=="chaotic"])

tab2=data.frame(TSlength=sims_sub$TSlength,NoiseLevel=sims_sub$NoiseLevel, Classification=sims_sub$Classification, LEtest=LEpoint)
err=aggregate(LEtest~NoiseLevel*TSlength*Classification, data=tab2, FUN=function(x) {length(which(x>0))/length(x)})
err$NoiseLevel2=ifelse(err$NoiseLevel==0, 0.01, err$NoiseLevel)

ggplot(err, aes(x=factor(NoiseLevel2), y=factor(TSlength), fill=LEtest)) +
  facet_grid(.~Classification) +
  geom_tile() + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  geom_text(aes(label=round(LEtest,2)), color="black", size=3) +
  labs(y="Time Series Length", x="Noise Level", fill="Prop chaotic") +
  scale_fill_distiller(palette = "Blues", direction = 1, limits=c(0,1))

# asymptotic SE method ####

LEasympresultsall=map2(modelresultsbest,sims_results_sub$jacobians,LEasymp)
#save(LEasympresultsall, file = "data/LEasympresultsall.RData")

#get CI
LEasall=map_dbl(LEasympresultsall, function(x) x$LElower)

#FNR
length(which(LEasall[sims_sub$Classification=="chaotic"]<0))/length(LEasall[sims_sub$Classification=="chaotic"])
#FPR
length(which(LEasall[sims_sub$Classification!="chaotic"]>0))/length(LEasall[sims_sub$Classification!="chaotic"])

tab2=data.frame(TSlength=sims_sub$TSlength,NoiseLevel=sims_sub$NoiseLevel, Classification=sims_sub$Classification, LEtest=LEasall)
fpr=aggregate(LEtest~NoiseLevel*TSlength*Classification, data=tab2, FUN=function(x) {length(which(x>0))/length(x)})
fpr$NoiseLevel2=ifelse(fpr$NoiseLevel==0, 0.01, fpr$NoiseLevel)

ggplot(fpr, aes(x=factor(NoiseLevel2), y=factor(TSlength), fill=LEtest)) +
  facet_grid(.~Classification) +
  geom_tile() + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  geom_text(aes(label=round(LEtest,2)), color="black", size=3) +
  labs(y="Time Series Length", x="Noise Level", fill="Prop chaotic") +
  scale_fill_distiller(palette = "Blues", direction = 1, limits=c(0,1))

#importance sampling method ####

#testing
# testsd=smap_model(data = sims_sub$data[[1]],y = "Value", ylog=F,pgr = "none",
#                   Efix = 2, taufix = sims_sub$taubest[1])
# st=Sys.time()
# test2=LEisamp(modelresultsbest[[2]], nreps = 1000)
# en=Sys.time(); en-st
# 
# test2$JLE[which.min(abs(test2$Wcum-0.025))]
# test2$JLE[which.min(abs(test2$Wcum-0.975))]

#############
# #batch run 
# 9 hours for 6 sim reps
# #free up some memory
# rm(sims_d, sims_results)
# plan(multisession, workers = 5)
# st=Sys.time()
# LEisampresults4=future_map(modelresultsbest,LEisamp,nreps=1000)
# en=Sys.time(); en-st
# save(LEisampresults4, file = "data/LEisampresults4.RData")
# 
# st2=Sys.time()
# LEbootresults4=future_map(modelresultsbest,LEboot,nreps=1000)
# en2=Sys.time(); en2-st2; en2-st
# plan(sequential)
# save(LEbootresults4, file = "data/LEbootresults4.RData")
# 
# load("data/LEbootresults.RData")
# load("data/LEbootresults2.RData")
# load("data/LEbootresults3.RData")
# load("data/LEisampresults.RData")
# load("data/LEisampresults2.RData")
# load("data/LEisampresults3.RData")
# 
# LEbootresults_all=c(LEbootresults,LEbootresults2,LEbootresults3,LEbootresults4)
# LEisampresults_all=c(LEisampresults,LEisampresults2,LEisampresults3,LEisampresults4)
# save(LEbootresults_all, file = "data/LEbootresults_all.RData")
# save(LEisampresults_all, file = "data/LEisampresults_all.RData")
#############

plan(multisession, workers = 5)
st=Sys.time()
LEisampresults_all=future_map(modelresultsbest,LEisamp,nreps=1000)
en=Sys.time(); en-st
plan(sequential)
#save(LEisampresults, file = "data/LEisampresults.RData")

#get CI 
LEis=map_dbl(LEisampresults_all, function(x) x$JLE[which.min(abs(x$Wcum-0.05))])

#FPR
length(which(LEis[sims_sub$Classification!="chaotic"]>0))/length(LEis[sims_sub$Classification!="chaotic"])
#FNR
length(which(LEis[sims_sub$Classification=="chaotic"]<0))/length(LEis[sims_sub$Classification=="chaotic"])

tab2=data.frame(TSlength=sims_sub$TSlength,NoiseLevel=sims_sub$NoiseLevel, Classification=sims_sub$Classification, LEtest=LEis)
err=aggregate(LEtest~NoiseLevel*TSlength*Classification, data=tab2, FUN=function(x) {length(which(x>0))/length(x)})
err$NoiseLevel2=ifelse(err$NoiseLevel==0, 0.01, err$NoiseLevel)

ggplot(err, aes(x=factor(NoiseLevel2), y=factor(TSlength), fill=LEtest)) +
  facet_grid(.~Classification) +
  geom_tile() + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  geom_text(aes(label=round(LEtest,2)), color="black", size=3) +
  labs(y="Time Series Length", x="Noise Level", fill="Prop chaotic") +
  scale_fill_distiller(palette = "Blues", direction = 1, limits=c(0,1))

#residual bootstrap method ####
#get residuals for response variable (gr), resample, add to predicted
#refit model with same hyperpars but new response var and get new theta, same predictors
#get gle

#plan(multisession, workers = 5)
st=Sys.time()
LEbootresults_all=future_map(modelresultsbest,LEboot,nreps=1000)
en=Sys.time(); en-st
plan(sequential)
#save(LEbootresults, file = "data/LEbootresults.RData")

#get CI 
LEbt=map_dbl(LEbootresults_all,quantile,0.05)

#FPR
length(which(LEbt[sims_sub$Classification!="chaotic"]>0))/length(LEbt[sims_sub$Classification!="chaotic"])
#FNR
length(which(LEbt[sims_sub$Classification=="chaotic"]<0))/length(LEbt[sims_sub$Classification=="chaotic"])

tab2=data.frame(TSlength=sims_sub$TSlength,NoiseLevel=sims_sub$NoiseLevel, Classification=sims_sub$Classification, LEtest=LEbt)
err=aggregate(LEtest~NoiseLevel*TSlength*Classification, data=tab2, FUN=function(x) {length(which(x>0))/length(x)})
err$NoiseLevel2=ifelse(err$NoiseLevel==0, 0.01, err$NoiseLevel)

ggplot(err, aes(x=factor(NoiseLevel2), y=factor(TSlength), fill=LEtest)) +
  facet_grid(.~Classification) +
  geom_tile() + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  geom_text(aes(label=round(LEtest,2)), color="black", size=3) +
  labs(y="Time Series Length", x="Noise Level", fill="Prop chaotic") +
  scale_fill_distiller(palette = "Blues", direction = 1, limits=c(0,1))
