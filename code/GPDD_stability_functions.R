# Functions for GPDD stability analysis

#get max E, theta, R2, and s-map coeffs for each series

#true R2 function
getR2=function(resultsdf, obse="obs", prede="pred") {
  d=na.omit(select(resultsdf, obse, prede))
  R2=1-sum((d[,obse]-d[,prede])^2)/sum((d[,obse]-mean(d[,obse]))^2)
}

besthyper=function(data=NULL, Efix=NULL, taufix=NULL, y, ylog, pgr=NULL, returntable=F) { 
  #get best E, tau, theta
  #y is the time series
  if(ylog==T & pgr=="fd") stop("must use ylog=F for pgr='fd'")
  #if providing dataframe containing y, give y as character (column name)
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
  #set max E
  #Emax=min(12,round(sqrt(serlen),0))
  #Emax=min(12,floor(sqrt(serlen))) #round down
  
  if(!is.null(taufix)) {
    tautry=taufix #set tau to fixed value
  } else {
    #tautry=1:Emax
    tautry=1:6 #24
  }
  if(!is.null(Efix)) {
    Etry=Efix #set E to fixed value
  } else {
    #Etry=1:Emax
    Etry=1:6 #12
  }
  
  #thetatest = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
  thetatest = c(0, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
  
  # model_selection=expand.grid(Etry, tautry)
  # colnames(model_selection)=c("E","tau")
  
  model_selection=expand.grid(Etry, tautry, thetatest)
  colnames(model_selection)=c("E","tau","theta")
  #determine sort order
  model_selection=arrange(model_selection, tau, theta, E)
  
  #filter candidate models 
  #model_selection=filter(model_selection, E*tau<=Emax)
  if(length(Etry)!=1 & length(tautry)!=1) {
    #if E or tau are both free
    #model_selection=filter(model_selection, E*(E+tau)<serlen & E*tau/serlen<0.2)
    model_selection=filter(model_selection, E^2<=serlen & E*tau/serlen<=0.2)
  } else if(length(Etry)!=1 | length(tautry)!=1) {
    #if one of E or tau are free, the other fixed
    #model_selection=filter(model_selection, E*(E+tau)<serlen)
    model_selection=filter(model_selection, E^2<=serlen)
  } else {
    #if both fixed, make sure E isn't too large
    #while(model_selection$E*(model_selection$E+model_selection$tau)>serlen) {
    while(model_selection$E^2>serlen) {
        model_selection$E=model_selection$E-1
    }
  }
  model_selection$rmse=NA
  model_selection$npred=NA
  
  #generate lags
  #ser_lags=make_block(ser, max_lag = Emax+1, tau=1)[,-1]
  ser_lags=make_block(ser, max_lag = max(model_selection$E*model_selection$tau)+1, tau=1)[,-1]

  #ser_lags=na.omit(ser_lags) #standardize targets
  #get best E and tau using simplex
  for(i in 1:nrow(model_selection)) {
    #if using growth rate as response, substitute it
    if(pgr=="fd") {
      ser_lags$col1=ser-lag(ser,model_selection$tau[i])
    } 
    if(pgr=="gr") {
      ser_lags$col1=ser_log-lag(ser_log,model_selection$tau[i])
    } 
    # simptemp=block_lnlp(ser_lags, tp = 0, method = "simplex", 
    #                     columns = paste0("col1_",model_selection$tau[i]*(1:model_selection$E[i])), 
    #                     target_column = "col1",
    #                     first_column_time = F, silent = T, stats_only = F)
    simptemp=block_lnlp(ser_lags, tp = 0, method = "s-map", 
                        columns = paste0("col1_",model_selection$tau[i]*(1:model_selection$E[i])), 
                        target_column = "col1", theta=model_selection$theta[i],
                        first_column_time = F, silent = T, stats_only = F)
    model_selection$rmse[i]=simptemp$rmse
    model_selection$npred[i]=simptemp$num_pred
    resultsdf=simptemp$model_output[[1]]
    resultsdf$Obs_abund=ser_or
    #get error based on abundance
    if(pgr=="none" & ylog==F) {
      resultsdf$Pred_abund=resultsdf$pred
    }  
    if(pgr=="none" & ylog==T) {
      resultsdf$Pred_abund=exp(resultsdf$pred)
    }
    if(pgr=="fd") {
      resultsdf$Pred_abund=lag(ser_or,model_selection$tau[i])+resultsdf$pred
      #model_selection$rmse[i]=sqrt(sum((ser_or-Pred_abund)^2, na.rm=T)/model_selection$npred[i])
    }
    if(pgr=="gr") {
      resultsdf$Obs_abund=ser_or
      resultsdf$Pred_abund=lag(ser_or,model_selection$tau[i])*exp(resultsdf$pred)
      #model_selection$rmse[i]=sqrt(sum((ser_or-Pred_abund)^2, na.rm=T)/model_selection$npred[i])
    }
    model_selection$rmse[i]=-getR2(resultsdf, obse="Obs_abund", prede="Pred_abund")
  }
  #calculate error accounting for df
  #model_selection$error=with(model_selection, ifelse(npred-E^2<0, NA, npred*rmse^2/(npred-E^2)))
  model_selection$error=round(model_selection$rmse,2)
  bestE=model_selection$E[which.min(model_selection$error)]
  bestTau=model_selection$tau[which.min(model_selection$error)]
  bestTheta=model_selection$theta[which.min(model_selection$error)]
  
  # #get best theta
  # ser_lags=make_block(ser, max_lag = bestE+1, tau=bestTau)[,-1]
  # #if using growth rate as response, substitute it
  # if(pgr=="fd") {
  #   ser_lags$col1=ser-lag(ser,bestTau)
  # } 
  # if(pgr=="gr") {
  #   ser_lags$col1=ser_log-lag(ser_log,bestTau)
  # } 
  # smap_results=block_lnlp(ser_lags, tp = 0, method = "s-map", 
  #                         columns = paste0("col1_",bestTau*(1:bestE)),
  #                         target_column = "col1", theta=thetatest,
  #                         first_column_time = F, silent = T)
  # bestTheta=smap_results$theta[which.min(smap_results$rmse)]
  # 
  if(returntable) {
    return(model_selection)
  } else {
    return(data.frame(bestE=bestE, bestTau=bestTau, bestTheta=bestTheta, Emax=max(model_selection$E), taumax=max(model_selection$tau)))
  }
}

smap_model=function(data, hypars=NULL, y, ylog, pgr="none", Efix=NULL, taufix=NULL) {
  #y col name for untransformed data
  #ylog indicates whether to y log transform y
  #to use different response, set pgr to "gr" (growth rate) or "fd" (first difference)
  #if pgr="fd", must use ylog=F
  #get hyper parameters
  if(!is.null(hypars)) {
    hyperpars=hypars
  } else {
    hyperpars=besthyper(data=data, y=y, ylog=ylog, pgr=pgr, Efix=Efix, taufix=taufix)
  }
  ser_or=data[,y]
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
  return(list(modelstats=modelstats, resultsdf=resultsdf, form=form))
}

#select among multiple models
smap_model_options=function(data, hypars=NULL, Efix=NULL, taufix=NULL, y, model) {
  if(model==1) {
    modelresults=smap_model(data, hypars=hypars, y=y, ylog=F, Efix=Efix, taufix=taufix)
  }
  if(model==2) {
    modelresults=smap_model(data, hypars=hypars, y=y, ylog=T, Efix=Efix, taufix=taufix)
  }
  if(model==3) {
    modelresults=smap_model(data, hypars=hypars, y=y, pgr="fd", ylog=F, Efix=Efix, taufix=taufix)
  }
  if(model==4) {
    modelresults=smap_model(data, hypars=hypars, y=y, pgr="gr", ylog=F, Efix=Efix, taufix=taufix)
  }
  if(model==5) {
    modelresults=smap_model(data, hypars=hypars, y=y, pgr="gr", ylog=T, Efix=Efix, taufix=taufix)
  }
  return(modelresults)
}

#construct Jacobians
getJacobians=function(modelresults){
  #this will vary depending on the model, data transform
  form=modelresults$form
  
  ndim=modelresults$modelstats$E #dimension
  tau=modelresults$modelstats$tau #tau
  len=nrow(modelresults$resultsdf)
  coefs=modelresults$resultsdf[,paste0("c_",1:ndim)]
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

#get local and global LEs
getStability=function(modelresults, jacobians) {
  tau=modelresults$modelstats$tau
  
  len=dim(jacobians)[3]
  ndim=dim(jacobians)[1]
  
  if(ndim==1) {
    lle=log(abs(jacobians[1,1,]))/tau
  } else {
    lle=numeric(length = len)
    for(i in 1:len) {
      Jac=jacobians[,,i]
      if(any(is.na(Jac))) {
        lle[i]=NA
      } else {
        lle[i]=log(max(abs(eigen(Jac, only.values = T)$values)))/tau
      }
    }
  }
  
  #proportion positive
  lle_pp=length(which(lle>0))/length(which(!is.na(lle)))
  
  #naive avg of local LEs
  lle_avg=mean(lle, na.rm=T)
  
  #global le
  runlen=rle(!is.na(jacobians[1,1,]))
  serlen=max(runlen$lengths[runlen$values==TRUE])
  if(serlen<(len/2.5) & serlen<25) {
    #do not compute if the longest run of non-missing values is less than 40% of ts length or less than 25
    gle=NA 
  } else {
    if(ndim==1) {
      #average over all timepoints
      gle=mean(log(abs(jacobians[1,1,])), na.rm=T)/tau
    } else {
      #compute over longest run of non-missing values
      starti=sum(runlen$lengths[1:(which.max(runlen$lengths)-1)])+1
      endi=sum(runlen$lengths[1:(which.max(runlen$lengths))])
      LEtemp=numeric(tau)
      for(i in 1:tau){
        indices=seq(from=starti+i-1, to=endi, by=tau)
        Jacs=jacobians[,,indices]
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
        LEtemp[i]=1/nk*log(max(abs(diag(Rcum))))/tau
      }
      gle=mean(LEtemp)
    }
  }
  
  return(list(lle=lle, lle_pp=lle_pp, lle_avg=lle_avg, gle=gle))
}

#1d LE from best model
LE1d=function(data, bestmodel, taufix=NULL, y) {
  modelresults=smap_model_options(data=data, y=y, model=bestmodel, Efix=1, taufix=taufix)
  jacobians=getJacobians(modelresults)
  stability=getStability(modelresults, jacobians)
  return(list(modelresults=modelresults, jacobians=jacobians, stability=stability))
}

# Regression-based LE (Rosenstein method) 
regLE=function(data, modelresults, y) {
  #recompute bestE with simplex
  ser=data[,y]
  bestE=modelresults$modelstats$E #dimension
  tau=modelresults$modelstats$tau #tau
  
  xi=make_block(ser, max_lag = bestE, tau=tau)[-(1:((bestE-1)*tau)),-1] #generate E lags
  if(bestE==1) {ni=length(xi)} else {ni=nrow(xi)}
  navals=which(!complete.cases(xi))
  rownames(xi)=NULL
  if(bestE==1) {xi[navals]=NA} else {xi[navals,]=NA}
  
  #distij=as.matrix(dist(xi))+10000*diag(ni) 
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
  # #drop initial distance
  # xr=matrix(ncol=2, nrow=smax, data = c(rep(1,smax), 1:smax))
  # avg_dist=avg_dist[2:(smax+1)]
  bb=solve(t(xr)%*%xr)%*%t(xr)%*%avg_dist
  verr=t(avg_dist-xr%*%bb)%*%(avg_dist-xr%*%bb)/(smax-1)
  vcoef=verr[1,1]*solve(t(xr)%*%xr)
  # reg=lm(avg_dist~c(0:smax))
  # summary(reg)
  # plot(0:smax, avg_dist)
  
  LEreg=bb[2]
  LEreg_se=sqrt(vcoef[2,2])
  
  return(data.frame(LEreg=LEreg, LEreg_se=LEreg_se))
}

#gets LE by averaging segments 
LEshift=function(modelresults, jacobians, samplinginterval="monthly") {
  tau=modelresults$modelstats$tau

  len=dim(jacobians)[3]
  ndim=dim(jacobians)[1]
  
  runlen=rle(!is.na(jacobians[1,1,]))
  serlen=max(runlen$lengths[runlen$values==TRUE])
  
  Tminus=3:6
  LEseg=data.frame(SegLen=(serlen-max(Tminus)):(serlen-min(Tminus)), le_mean=NA, le_sd=NA, le_ci=NA, le_n=NA, le_var=NA, le_var_mo=NA) %>% 
    filter(SegLen>0)
  
  for(i in 1:nrow(LEseg)) {
    SegLen=LEseg$SegLen[i]
    LEtemp=numeric(len-SegLen+1)
    for(j in 1:length(LEtemp)) {
      Jacs1=jacobians[,,j:(j+SegLen-1)]
      LEtemp2=numeric(tau)
      for(a in 1:tau) {
        indices=seq(from=1+a-1, to=ifelse(ndim==1,length(Jacs1),dim(Jacs1)[3]), by=tau)
        if(ndim==1) {
          Jacs=Jacs1[indices]
        } else {
          Jacs=Jacs1[,,indices]
        }
        if(any(is.na(Jacs))) {
          LEtemp2[a]=NA
        } else if(ndim==1) {
          LEtemp2[a]=mean(log(abs(Jacs)), na.rm=T)/tau
        } else {
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
          LEtemp2[a]=1/nk*log(max(abs(diag(Rcum))))/tau
        }
      }
      LEtemp[j]=mean(LEtemp2, na.rm=T)
    }
    LEseg$le_n[i]=length(which(!is.na(LEtemp)))
    LEseg$le_mean[i]=mean(LEtemp, na.rm=T)
    LEseg$le_sd[i]=sd(LEtemp, na.rm=T)
    LEseg$le_ci[i]=LEseg$le_sd[i]/sqrt(LEseg$le_n[i])*qt(p=0.95, df=LEseg$le_n[i]-1)
    LEseg$le_var[i]=var(LEtemp, na.rm=T)
    LEseg$le_var_mo[i]=var(LEtemp/timescale_mo(samplinginterval, 1), na.rm=T)
  }
  
  minmean=min(LEseg$le_mean) #not necessarily the mean with lowest CI
  minci=min(LEseg$le_mean-LEseg$le_ci)
  varmin=min(LEseg$le_var)
  varmin_mo=min(LEseg$le_var_mo)

  return(list(LEseg=LEseg, minmean=minmean, minci=minci, varmin=varmin, varmin_mo=varmin_mo))
}

#converts timesteps to months
timescale_mo=function(SamplingInterval, y) {
  ifelse(SamplingInterval=="annual", y*12,
         ifelse(SamplingInterval=="seasonal", y*6,
                ifelse(SamplingInterval=="bimonthly", y*2,
                       ifelse(SamplingInterval=="monthly", y,
                              ifelse(SamplingInterval=="4-week", y/1.07,
                                     ifelse(SamplingInterval=="weekly", y/4.286,NA))))))
  
}

SiblymodelLE=function(data, y) {
  ser=data[,y]
  ser_log=log(ser)
  p1=lag(ser_log)
  p2=lag(ser_log)^2
  mod=lm(ser_log~p1+p2)
  c1=coefficients(mod)[2]
  c2=coefficients(mod)[3]
  LE=mean(log(abs(c1+2*c2*ser_log)), na.rm=T)
}