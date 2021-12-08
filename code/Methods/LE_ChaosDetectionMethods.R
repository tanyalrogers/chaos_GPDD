# This script contains several functions that are used to:
# - select the best E, tau, and theta
# - get Lyaponov exponents (LEs) using the Jacobian method (B2 in paper)
# - get LEs using the direct method (B1 in paper)

#### Jacobian LE Method (B2) ####

#R2 function
getR2=function(resultsdf, obse="obs", prede="pred") {
  d=na.omit(select(resultsdf, obse, prede))
  R2=1-sum((d[,obse]-d[,prede])^2)/sum((d[,obse]-mean(d[,obse]))^2)
}

#Get best hyperparameters (E, tau, theta)
besthyper=function(
  data=NULL, #a data frame containing the time series
  Efix=NULL, #option to fix E to a particular value
  taufix=NULL, #option to fix tau to a particular value
  y, #either a time series vector, or the name of the time series variable in 'data'
  ylog, #whether or not to log the data (T/F)
  pgr="none", #whether response variable should be untransformed "none", first difference "fd", or growth rate "gr"
  returntable=F, #return full model selection table rather than just best model
  yboot=NULL #bootstrapped response variable
  ) { 
  
  #if using first difference as response, ylog must be set to F
  if(ylog==T & pgr=="fd") stop("must use ylog=F for pgr='fd'")
  
  #if providing data frame containing y, give y as character (column name)
  if(!is.null(data)) {
    ser_or=data[,y]
  } else {
    ser_or=y
  }
  
  #log transform
  if(ylog==T) ser=log(ser_or) else ser=ser_or
  if(pgr=="gr") ser_log=log(ser_or)
  
  #get effective ts length (longest string of non-NA values)
  len=rle(!is.na(ser))
  serlen=max(len$lengths[len$values==TRUE])
  
  #range of E and tau to try
  #max for both is currently hard coded as 6, but could be changed/added as an argument 
  if(!is.null(taufix)) {
    tautry=taufix #set tau to fixed value
  } else {
    tautry=1:6
  }
  if(!is.null(Efix)) {
    Etry=Efix #set E to fixed value
  } else {
    Etry=1:6
  }
  
  #range of theta to try
  thetatest = c(0, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
  
  #create table of candidate models
  model_selection=expand.grid(Etry, tautry, thetatest)
  colnames(model_selection)=c("E","tau","theta")
  #determine sort order
  model_selection=arrange(model_selection, tau, theta, E)
  
  #filter candidate models to reasonable E and tau combinations given effective ts length
  if(length(Etry)!=1 & length(tautry)!=1) {
    #if E or tau are both free
    model_selection=filter(model_selection, E^2<=serlen & E*tau/serlen<=0.2)
  } else if(length(Etry)!=1 | length(tautry)!=1) {
    #if one of E or tau are free, the other fixed
    model_selection=filter(model_selection, E^2<=serlen)
  } else {
    #if both fixed, make sure E isn't too large
    while(unique(model_selection$E)^2>serlen) {
        model_selection$E=model_selection$E-1
    }
  }
  model_selection$rmse=NA
  model_selection$npred=NA
  model_selection$error1=NA
  
  #generate lags
  ser_lags=make_block(ser, max_lag = max(model_selection$E*model_selection$tau)+1, tau=1)[,-1]
  #ser_lags=na.omit(ser_lags) #if you wish to standardize targets across models
  
  #run s-map for all parameter sets
  for(i in 1:nrow(model_selection)) {
    #if using first difference or growth rate as response, substitute it
    if(pgr=="fd") {
      ser_lags$col1=ser-lag(ser,model_selection$tau[i])
    } 
    if(pgr=="gr") {
      ser_lags$col1=ser_log-lag(ser_log,model_selection$tau[i])
    } 
    if(!is.null(yboot)) { #if bootstrap, swap response variable
      ser_lags$col1=yboot
    }
    simptemp=block_lnlp(ser_lags, tp = 0, method = "s-map", 
                        columns = paste0("col1_",model_selection$tau[i]*(1:model_selection$E[i])), 
                        target_column = "col1", theta=model_selection$theta[i],
                        first_column_time = F, silent = T, stats_only = F)
    model_selection$rmse[i]=simptemp$rmse
    model_selection$npred[i]=simptemp$num_pred
    resultsdf=simptemp$model_output[[1]]
    resultsdf$Obs_abund=ser_or
    #get error based on abundance
    #requires backtransformation if using first difference or growth rate
    if(pgr=="none" & ylog==F) {
      resultsdf$Pred_abund=resultsdf$pred
    }  
    if(pgr=="none" & ylog==T) {
      resultsdf$Pred_abund=exp(resultsdf$pred)
    }
    if(pgr=="fd") {
      resultsdf$Pred_abund=lag(ser_or,model_selection$tau[i])+resultsdf$pred
    }
    if(pgr=="gr") {
      resultsdf$Pred_abund=lag(ser_or,model_selection$tau[i])*exp(resultsdf$pred)
    }
    model_selection$error1[i]=-getR2(resultsdf, obse="Obs_abund", prede="Pred_abund")
  }
  
  #calculate error accounting for df
  model_selection$error=round(model_selection$error1,2)
  bestE=model_selection$E[which.min(model_selection$error)]
  bestTau=model_selection$tau[which.min(model_selection$error)]
  bestTheta=model_selection$theta[which.min(model_selection$error)]
  
  #output
  if(returntable) {
    return(model_selection)
  } else {
    return(data.frame(bestE=bestE, bestTau=bestTau, bestTheta=bestTheta, Emax=max(model_selection$E), taumax=max(model_selection$tau)))
  }
}

#Get s-map model output
#This uses the above function to get the best hyperparameters.
smap_model=function(
  data=NULL, #a data frame containing the time series
  hypars=NULL, #optional output of 'besthyper' containing pre-optimized hyperparameters 
  y, #either a time series vector, or the name of the time series variable in 'data'
  ylog, #whether or not to log the data (T/F)
  pgr="none", #whether response variable should be untransformed "none", first difference "fd", or growth rate "gr"
  Efix=NULL, #option to fix tau to a particular value
  taufix=NULL, #option to fix tau to a particular value
  yboot=NULL #bootstrapped response variable
  ) {

  #get hyper parameters if not already supplied
  #obviously, if you fix E and tau, this will return only a single model
  if(!is.null(hypars)) {
    hyperpars=hypars
  } else {
    hyperpars=besthyper(data=data, y=y, ylog=ylog, pgr=pgr, Efix=Efix, taufix=taufix, yboot=yboot)
  }
  
  #if providing data frame containing y, give y as character (column name)
  if(!is.null(data)) {
    ser_or=data[,y]
  } else {
    ser_or=y
  }
  
  #log transform
  if(ylog==T) ser=log(ser_or) else ser=ser_or
  if(pgr=="gr") ser_log=log(ser_or)
  
  #store smap output
  ser_lags=make_block(ser, max_lag = hyperpars$bestE+1, tau=hyperpars$bestTau)[,-1]
  if(pgr=="fd") {
    ser_lags$col1=ser-lag(ser,hyperpars$bestTau)
  } 
  if(pgr=="gr") {
    ser_lags$col1=ser_log-lag(ser_log,hyperpars$bestTau)
  } 
  if(!is.null(yboot)) { #if bootstrap, swap response variable
    ser_lags$col1=yboot
  }
  smap_results=block_lnlp(ser_lags, tp = 0, method = "s-map", 
                          columns = paste0("col1_",hyperpars$bestTau*(1:hyperpars$bestE)), 
                          target_column = "col1", theta=hyperpars$bestTheta,
                          first_column_time = F, silent = T,
                          stats_only = F, save_smap_coefficients = T)
  resultsdf=cbind(smap_results$model_output[[1]], smap_results$smap_coefficients[[1]])
  
  #untransformed abundance
  resultsdf$Obs_abund=ser_or
  #get untranformed abundance predictions based on model type
  if(pgr=="none" & ylog==F) {
    resultsdf$Pred_abund=resultsdf$pred
    form="ut-ut"
  }  
  if(pgr=="none" & ylog==T) {
    resultsdf$Pred_abund=exp(resultsdf$pred)
    form="log-log"
  }
  if(pgr=="fd") {
    resultsdf$Pred_abund=lag(resultsdf$Obs_abund,hyperpars$bestTau)+resultsdf$pred
    form="fd-ut"
  }
  if(pgr=="gr") {
    resultsdf$Pred_abund=lag(resultsdf$Obs_abund,hyperpars$bestTau)*exp(resultsdf$pred)
    #ylog is not used to get predictions, but used for jacobian construction
    if(ylog==T) form="gr-log"
    if(ylog==F) form="gr-ut"
  }
  
  #R2 for model
  modelR2=getR2(resultsdf) 
  #R2 for untransformed abundance
  modelR2_abund=getR2(resultsdf, obse="Obs_abund", prede="Pred_abund")
  modelstats=data.frame(E=hyperpars$bestE, tau=hyperpars$bestTau, theta=hyperpars$bestTheta, 
                        Emax=hyperpars$Emax, taumax=hyperpars$taumax, num_pred=smap_results$num_pred,
                        R2model=modelR2, R2abund=modelR2_abund, rho=smap_results$rho)
  return(list(modelstats=modelstats, resultsdf=resultsdf, form=form, smap_results=smap_results, ser_lags=ser_lags))
}

#Calls smap_model with several pre-specified model forms
#Note that models 1&3 and 2&5 produce identical (backtransformed) output
smap_model_options=function(data, hypars=NULL, Efix=NULL, taufix=NULL, y, model, yboot=NULL) {
  if(model==1) { #"ut-ut", response and predictors both untransformed
    modelresults=smap_model(data, hypars=hypars, y=y, ylog=F, Efix=Efix, taufix=taufix, yboot=yboot)
  }
  if(model==2) { #"log-log", response and predictors both log transformed
    modelresults=smap_model(data, hypars=hypars, y=y, ylog=T, Efix=Efix, taufix=taufix, yboot=yboot)
  }
  if(model==3) { #"fd-ut", response = first difference, predictors untransformed
    modelresults=smap_model(data, hypars=hypars, y=y, pgr="fd", ylog=F, Efix=Efix, taufix=taufix, yboot=yboot)
  }
  if(model==4) { #"gr-ut", response = growth rate, predictors untransformed
    modelresults=smap_model(data, hypars=hypars, y=y, pgr="gr", ylog=F, Efix=Efix, taufix=taufix, yboot=yboot)
  }
  if(model==5) { #"gr-log", response = growth rate, predictors log transformed
    modelresults=smap_model(data, hypars=hypars, y=y, pgr="gr", ylog=T, Efix=Efix, taufix=taufix, yboot=yboot)
  }
  return(modelresults)
}

#Construct Jacobian matrices from smap coefficients
getJacobians=function(
  modelresults #output from smap_model or smap_model_options
  ){
  
  #jacobians are derived in terms of untransformed abudance, so will vary depending on the model and data transform.
  #the LE is (in theory) not affected by the particular form of the model or Jacobian, but this ensures that the local LEs are correct.
  form=modelresults$form
  
  ndim=modelresults$modelstats$E #dimension
  tau=modelresults$modelstats$tau #tau
  len=nrow(modelresults$resultsdf)
  coefs=modelresults$resultsdf[,paste0("c_",1:ndim)]
  
  #jacobians are stored in a 3d array
  jacobians=array(dim = c(ndim,ndim,len))
  
  if(form=="ut-ut") {
    if(ndim==1) {
      jacobians[,,]=coefs
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1=matrix(nrow = ndim, ncol = ndim, 0)
          J1[1,]=as.numeric(coefs[i,])
          J1[2:ndim,1:(ndim-1)]=diag(ndim-1)
          jacobians[,,i]=J1
        }
      }
    }
  }
  
  if(form=="fd-ut") {
    if(ndim==1) {
      jacobians[,,]=coefs+1
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1=matrix(nrow = ndim, ncol = ndim, 0)
          J1[1,]=as.numeric(coefs[i,])+c(1,rep(0,ndim-1))
          J1[2:ndim,1:(ndim-1)]=diag(ndim-1)
          jacobians[,,i]=J1
        }
      }
    }
  }
  
  if(form=="log-log") {
    x_obs=modelresults$resultsdf$Obs_abund
    r_x=exp(modelresults$resultsdf$pred)
    if(ndim==1) {
      jacobians[,,]=(coefs)*r_x/lag(x_obs, tau)
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1<-J2<-J3<-matrix(nrow = ndim, ncol = ndim, 0)
          diag(J1)<-c(r_x[i], x_obs[i-(1:(ndim-1))*tau])
          #c(r_x[i], lag(x_obs)[i], lag(x_obs,2)[i], lag(x_obs,3)[i], lag(x_obs,4)[i])
          J2[1,]=as.numeric(coefs[i,])
          J2[2:ndim,1:(ndim-1)]=diag(ndim-1)
          diag(J3)<-1/x_obs[i-(1:(ndim))*tau] 
          #c(1/(lag(x_obs)[i]), 1/(lag(x_obs,2)[i]), 1/(lag(x_obs,3)[i]), 1/(lag(x_obs,4)[i]), 1/(lag(x_obs,5)[i]))
          jacobians[,,i]=J1%*%J2%*%J3
        }
      }
    }
  }
  
  if(form=="gr-ut") {
    x_obs=modelresults$resultsdf$Obs_abund
    r_x=exp(modelresults$resultsdf$pred)
    if(ndim==1) {
      jacobians[,,]=(coefs*lag(x_obs,tau)+1)*r_x
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1<-J2<-matrix(nrow = ndim, ncol = ndim, 0)
          diag(J1)<-c(r_x[i], rep(1,ndim-1))
          #c(r_x[i], 1, 1, 1)
          J2[1,]=as.numeric(coefs[i,])*x_obs[i-tau]+c(1,rep(0,ndim-1))
          J2[2:ndim,1:(ndim-1)]=diag(ndim-1)
          jacobians[,,i]=J1%*%J2
        }
      }
    }
  }
  
  if(form=="gr-log") {
    x_obs=modelresults$resultsdf$Obs_abund
    r_x=exp(modelresults$resultsdf$pred)
    if(ndim==1) {
      jacobians[,,]=(coefs+1)*r_x
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1<-J2<-J3<-matrix(nrow = ndim, ncol = ndim, 0)
          diag(J1)<-c(r_x[i]*x_obs[i-tau], x_obs[i-(1:(ndim-1))*tau])
          #c(r_x[i]*lag(x_obs)[i], lag(x_obs)[i], lag(x_obs,2)[i], lag(x_obs,3)[i], lag(x_obs,4)[i])
          J2[1,]=as.numeric(coefs[i,])+c(1,rep(0,ndim-1))
          J2[2:ndim,1:(ndim-1)]=diag(ndim-1)
          diag(J3)<-1/x_obs[i-(1:(ndim))*tau] 
          #c(1/(lag(x_obs)[i]), 1/(lag(x_obs,2)[i]), 1/(lag(x_obs,3)[i]), 1/(lag(x_obs,4)[i]), 1/(lag(x_obs,5)[i]))
          jacobians[,,i]=J1%*%J2%*%J3
        }
      }
    }
  }
  
  return(jacobians)
}

#Get LE (lower confidence bound) from jacobian matrices
#This is done by averaging LEs computed over several long subsegments of the data
LEshift=function(
  modelresults, #output from smap_model or smap_model_options
  jacobians #output from getJacobians
  ) {
  
  #modelresults is just used to get tau
  tau=modelresults$modelstats$tau
  
  len=dim(jacobians)[3] #time series length
  ndim=dim(jacobians)[1] #E
  
  #remove leading NAs
  if(ndim==1) {
    jacobians2=jacobians[(ndim*tau+1):len]
    runlen=rle(!is.na(jacobians2))
    len2=length(jacobians2)
  } else {
    jacobians2=jacobians[,,(ndim*tau+1):len]
    runlen=rle(!is.na(jacobians2[1,1,]))
    len2=dim(jacobians2)[3]
  }
  
  #effective ts length (longest string of non-NA values)
  serlen=max(runlen$lengths[runlen$values==TRUE])
  
  Tminus=3:6 #segment lengths range in length from T-3 to T-6
  LEseg=data.frame(SegLen=(serlen-max(Tminus)):(serlen-min(Tminus)), le_mean=NA, le_sd=NA, le_ci=NA, le_n=NA) %>% 
    filter(SegLen>0)
  
  for(i in 1:nrow(LEseg)) {
    SegLen=LEseg$SegLen[i]
    LEtemp=numeric(len2-SegLen+1)
    for(j in 1:length(LEtemp)) {
      if(ndim==1) {
        Jacs1=jacobians2[j:(j+SegLen-1)]
      } else {
        Jacs1=jacobians2[,,j:(j+SegLen-1)]
      }
      LEtemp2=numeric(tau)
      if(any(is.na(Jacs1))) {
        LEtemp2=NA
      } else {
        for(a in 1:tau) {
          indices=seq(from=a, to=ifelse(ndim==1,length(Jacs1),dim(Jacs1)[3]), by=tau)
          if(ndim==1) {
            Jacs=Jacs1[indices]
            LEtemp2[a]=mean(log(abs(Jacs)), na.rm=T)/tau
          } else {
            Jacs=Jacs1[,,indices]
            nk=dim(Jacs)[3]
            Jk=Jacs[,,1]
            QR=qr(Jk)
            R=qr.R(QR)
            Q=qr.Q(QR)
            Rcum=R
            for(k in 2:nk) {
              Jk=Jacs[,,k]
              QR=qr(Jk%*%Q)
              R=qr.R(QR)
              Q=qr.Q(QR)
              Rcum=R%*%Rcum
            }
            LEtemp2[a]=log(max(abs(diag(Rcum))))/tau/nk
          }
        }
      }
      LEtemp[j]=mean(LEtemp2, na.rm=F) #if tau>1, will result in more than i+1 segments being averaged unless you set na.rm=F
    }
    LEseg$le_n[i]=length(which(!is.na(LEtemp)))
    LEseg$le_mean[i]=mean(LEtemp, na.rm=T)
    LEseg$le_sd[i]=sd(LEtemp, na.rm=T)
    LEseg$le_ci[i]=LEseg$le_sd[i]/sqrt(LEseg$le_n[i])*qt(p=0.95, df=LEseg$le_n[i]-1)
  }
  
  minmean=min(LEseg$le_mean) #this is not necessarily the mean with lowest CI
  minci=min(LEseg$le_mean-LEseg$le_ci) #lowest of the lower confidence bounds on LE
  
  return(list(LEseg=LEseg, minmean=minmean, minci=minci))
}

#Get LE from particular model with E fixed to 1
LE1d=function(data, model, taufix=NULL, y) {
  modelresults=smap_model_options(data=data, y=y, model=model, Efix=1, taufix=taufix)
  jacobians=getJacobians(modelresults)
  JLE=LEshift(modelresults, jacobians)
  return(list(modelresults=modelresults, jacobians=jacobians, JLE=JLE))
}

#Get LE from particular model with E and tau fixed
LEfix=function(data, model, Efix=NULL, taufix=NULL, y) {
  modelresults=smap_model_options(data=data, y=y, model=model, Efix=Efix, taufix=taufix)
  jacobians=getJacobians(modelresults)
  JLE=LEshift(modelresults, jacobians)
  return(list(modelresults=modelresults, JLE=JLE))
}

#### Direct LE Method (B1) ####

# Regression-based method of Rosenstein
regLE=function(
  data, #a data frame containing the time series
  modelresults, #output from smap_model or smap_model_options
  y #the name of the time series variable in 'data'
  ) {
  
  ser=data[,y]
  
  #modelresults is used to get E and tau
  #could be written to take E and tau as inputs
  bestE=modelresults$modelstats$E
  tau=modelresults$modelstats$tau 
  
  xi=make_block(ser, max_lag = bestE, tau=tau)[-(1:((bestE-1)*tau)),-1] #generate E lags
  if(bestE==1) {ni=length(xi)} else {ni=nrow(xi)}
  navals=which(!complete.cases(xi))
  rownames(xi)=NULL
  if(bestE==1) {xi[navals]=NA} else {xi[navals,]=NA}
  
  distij=as.matrix(dist(xi)) #distance matrix
  exclusion=1 #exclusion radius
  delta <- row(distij) - col(distij)
  distij[delta <= exclusion & delta >= -1*exclusion] <- distij[delta <= exclusion & delta >= -1*exclusion]+10000
  
  md=apply(distij,1,min,na.rm=T) #min distance vector
  md[navals]=NA
  minind=apply(distij,1,which.min) #nearest point index vector
  minind[navals]=NA
  smax=4 #max steps into future
  avg_dist=numeric(smax+1)
  md[md==0]<-NA #remove distances of 0 (otherwise fails when taking logs)
  avg_dist[1]=mean(log(md), na.rm=T) 
  for(s in 1:smax) {
    ind1=(1:ni)+s
    ind2=minind+s
    rows=which(ind1<=ni & ind2<=ni)
    if(bestE==1){
      ds=abs(xi[ind1[rows]]-xi[ind2[rows]])
    } else {
      ds=sqrt(rowSums((xi[ind1[rows],]-xi[ind2[rows],])^2))
    }
    ds[ds==0]<-NA #remove distances of 0 (otherwise fails when taking logs)
    avg_dist[s+1]=mean(log(ds), na.rm=T)
  }
  
  #regression
  xr=matrix(ncol=2, nrow=smax+1, data = c(rep(1,smax+1), 0:smax))
 
  # if you want to drop initial distance (0)
  # xr=matrix(ncol=2, nrow=smax, data = c(rep(1,smax), 1:smax))
  # avg_dist=avg_dist[2:(smax+1)]
  
  # plot(0:smax, avg_dist)
  
  bb=solve(t(xr)%*%xr)%*%t(xr)%*%avg_dist
  verr=t(avg_dist-xr%*%bb)%*%(avg_dist-xr%*%bb)/(smax-1)
  vcoef=verr[1,1]*solve(t(xr)%*%xr)
  
  #or if you prefer:
  # reg=lm(avg_dist~c(0:smax))
  # summary(reg)
  
  LEreg=bb[2]
  LEreg_se=sqrt(vcoef[2,2])
  
  return(data.frame(LEreg=LEreg, LEreg_se=LEreg_se))
}

#### Miscellaneous ####

#converts timesteps to months, given different sampling intervals
timescale_mo=function(SamplingInterval, y) {
  ifelse(SamplingInterval=="annual", y*12,
         ifelse(SamplingInterval=="seasonal", y*6,
                ifelse(SamplingInterval=="bimonthly", y*2,
                       ifelse(SamplingInterval=="monthly", y,
                              ifelse(SamplingInterval=="4-week", y/1.07,
                                     ifelse(SamplingInterval=="weekly", y/4.286,NA))))))
  
}
